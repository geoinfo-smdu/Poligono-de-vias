#!/usr/bin/env python
# coding: utf-8

import geopandas as gpd
import pandas as pd
import math
from shapely.geometry import Point, MultiLineString, LineString, MultiPoint
from shapely.ops import substring, nearest_points, split, linemerge
from rasterstats import zonal_stats
from datetime import datetime

dist = 12.00

df_s_s = gpd.read_file(f'gis/SIRGAS_SHP_distrito_polygon.shp')
df_l_s = gpd.read_file(f'gis/SIRGAS_SHP_logradouronbl.shp')

## Selecionando apenas um distrito para teste
# df_s_s = df_s_s[df_s_s.ds_codigo == '54']


for i in df_s_s.itertuples():
    df_s = df_s_s[df_s_s.ds_codigo == i.ds_codigo]

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print("****** início: ", dt_string)	
    print(f"Executando {list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].lower()}")
    
    df_pvias = gpd.read_file('resultado/poligono_de_vias.gpkg', layer=f"{list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].lower()}")

    df_s_sindex = df_s.sindex
    df_l_sindex = df_l_s.sindex

    # Get the bounding box coordinates of the Polygon as a list
    bounds = list(df_s.bounds.values[0])

    # Get the indices of the Points that are likely to be inside the bounding box of the given Polygon
    l_candidate_idx = list(df_l_sindex.intersection(bounds))
    l_candidates = df_l_s.loc[l_candidate_idx]

    # Logradouros na área de estudo
    df_l = gpd.clip(l_candidates, df_s)

    points = df_l.apply(lambda x: MultiPoint([x.geometry.interpolate(((i)*dist)+dist/2+(x.geometry.length%dist)/2) 
                                                for i in range(math.floor(x.geometry.length/dist))]), 
                            axis=1)

    df_points = gpd.GeoDataFrame(geometry=points)

    def ponto_maior_borda(g):
        linhas_de_borda = split(MultiPoint(g.coords[:]).convex_hull.boundary, MultiPoint(g.coords[:]))
        maior_borda = max(linhas_de_borda, key=(lambda x: x.length))
        ponto_maior_borda = Point(maior_borda.coords[0])
        return ponto_maior_borda

    df_pvias['geometry_ponto_inicial'] = gpd.GeoDataFrame(geometry=df_pvias.exterior.apply(lambda x: ponto_maior_borda(x)))

    def bordas(r):
        ponto_inicial = r.geometry.boundary.project(r.geometry_ponto_inicial, normalized=True)
        borda = r.geometry.boundary.difference(r.geometry_ponto_inicial.buffer(1))
        borda_merged = linemerge(borda) if borda.type != 'LineString' else borda
        try:
            bordas = MultiLineString([substring(borda_merged, 0, 0.5, normalized=True), substring(borda_merged, 0.5, 1, normalized=True)])
            return bordas
        except:
            print(borda_merged)

    df_pvias['geometry_bordas'] = df_pvias.apply(lambda x: bordas(x), axis=1)

    df_point_to_border = gpd.sjoin(df_points[~df_points.is_empty].reset_index(drop=True).explode(), df_pvias, op='within')

    def calc_distancia(p):
        try:
            l_dist = LineString([nearest_points(p.geometry, p.geometry_bordas[0])[1], 
                                nearest_points(p.geometry, p.geometry_bordas[1])[1]])
            if p.geometry.distance(l_dist) < 1:
                return l_dist
        except: 
            print(p)
        
    df_dists = gpd.GeoDataFrame(geometry=df_point_to_border.apply(lambda x: calc_distancia(x), axis=1))

    df_dists.reset_index(drop=True, inplace=True)

    df_dists['largura'] = df_dists.length

    df_dists.crs = 'EPSG:31983'

    df_largura_pontos = gpd.GeoDataFrame({'largura':df_dists.length}, geometry=df_dists.geometry.centroid)

    larguras = gpd.sjoin(df_pvias, df_largura_pontos, how='left', op='intersects')['largura'].groupby(level=0).agg(['min', 'max', 'mean', 'std', 'count'] )

    df_larguras = df_pvias.join(larguras).drop(columns=['geometry_ponto_inicial', 'geometry_bordas'])

    df_larguras.rename(columns={"min": "largura_minima",
                            "max": "largura_maxima",
                            "mean": "largura_media",
                            "std": "largura_desvio_padrao",
                            "count": "largura_contagem_amostras"},
                  inplace=True)

    raster_stats = zonal_stats(vectors=df_larguras['geometry'], 
                           raster='/media/fernando/0dcfe7e0-23bd-4899-9552-3722ef65582e/MDT_SPTOTAL/SP_MDT_las_5m.tif')

    df_larguras = df_larguras.join(pd.DataFrame(raster_stats))

    df_larguras.rename(columns={"min": "cota_minima",
                            "max": "cota_maxima",
                            "mean": "cota_media",
                            "count": "cota_contagem_amostras"},
                  inplace=True)

    df_larguras['declividade_maxima_estatistica_percentual'] =((df_larguras.cota_maxima - df_larguras.cota_minima) / 
              (df_larguras.area / df_larguras.largura_media)) * 100

    df_larguras.to_file(f"resultado/{list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].lower()}.gpkg", driver='GPKG')
    df_larguras.to_file('resultado/poligono_de_vias_com_larguras.gpkg', layer=f"{list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].lower()}", driver="GPKG", OVERWRITE='YES')