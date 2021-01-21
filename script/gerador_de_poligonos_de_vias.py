#!/usr/bin/env python
# coding: utf-8

# # Técnica de voronoi com pontos mais próximos do polígono
# 
# Nesse notebook vamos experimentar uma técnica nova para tentar corrigir o excesso de polígono gerado em cruzamentos de grandes vias

from shapely.geometry import MultiPoint, LineString, MultiLineString, Point
from shapely.ops import nearest_points
import numpy as np
from libpysal.cg.voronoi import voronoi, voronoi_frames

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt

from shapely.ops import unary_union, polygonize, cascaded_union

from shapely import affinity

from datetime import datetime

import os


# df_qf_s = gpd.read_file(f'gis/SIRGAS_SHP_quadraMDSF.shp')
df_qv_s = gpd.read_file(f'gis/SIRGAS_SHP_quadraviariaed_polygon.shp')
df_l_s = gpd.read_file(f'gis/SIRGAS_SHP_logradouronbl.shp')
# df_s = gpd.read_file(f'gis/SIRGAS_SHP_setorfiscal.shp')
# df_s = gpd.read_file(f'gis/SIRGAS_SHP_prefeitura_regional_polygon.shp')
df_s_s = gpd.read_file(f'gis/SIRGAS_SHP_distrito_polygon.shp')
# df_s = gpd.read_file(f'gis/teste-recorte-1.geojson')
df_represas = gpd.read_file(f'gis/SIRGAS_REPRESAS_NIVELMAX.shp')
df_massa_dagua = gpd.read_file(f'gis/SIRGAS_MASSADAGUA.shp')

tolerancia_do_cruzamento = 1.5 # era 0.6

# df_s_s = df_s_s[df_s_s['ds_nome'] == 'SE']

for i in df_s_s.itertuples():
    df_s = df_s_s[df_s_s.ds_codigo == i.ds_codigo]

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print("****** início: ", dt_string)	
    print(f"Executando {list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].lower()}")
    
    if os.path.exists(f"./resultado/{list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].lower()}"):
        continue

    # Data Frame dos canteiros
    df_canteiro = df_qv_s[df_qv_s['qe_tipo'] == 'Praca_Canteiro']

    # Guardando as quadras em uma variavel para uso futuro
    df_quadras = gpd.clip(df_qv_s[(df_qv_s['qe_tipo'] != 'Praca_Canteiro') & df_qv_s.is_valid], df_s)
    # df_quadras = gpd.clip(df_qf, df_s)

    # DAtaFRame das Quadras viárias que não são canteiros
    df_qv = df_qv_s[df_qv_s['qe_tipo'] != 'Praca_Canteiro']
    # df_qv = df_qf

    # REcorte dos canteiros na área de estudo
    df_ct = gpd.overlay(df_s, df_canteiro, how='intersection')

    # Removendo área de represa
    df_s = gpd.overlay(df_s, df_represas, how='difference')

    # Removendoo Massa dagua
    df_s = gpd.overlay(df_s, df_massa_dagua, how='difference')

    df_s_sindex = df_s.sindex
    df_l_sindex = df_l_s.sindex

    # Get the bounding box coordinates of the Polygon as a list
    bounds = list(df_s.bounds.values[0])

    # Get the indices of the Points that are likely to be inside the bounding box of the given Polygon
    l_candidate_idx = list(df_l_sindex.intersection(bounds))
    l_candidates = df_l_s.loc[l_candidate_idx]

    # Logradouros na área de estudo
    df_l = gpd.clip(l_candidates, df_s)

    df_qv_sindex = df_qv.sindex

    # Get the bounding box coordinates of the Polygon as a list
    bounds = list(df_s.bounds.values[0])

    # Get the indices of the Points that are likely to be inside the bounding box of the given Polygon
    qv_candidate_idx = list(df_qv_sindex.intersection(bounds))
    qv_candidates = df_qv.iloc[qv_candidate_idx]

    # Poligono geral de vias da área de estudo
    # Este que será recortado adequadamente
    df_vias = gpd.overlay(df_s, qv_candidates, how='difference')

    # Recortando apenas as vias que estão dentro do polígono de vias
    df_l_cropped = gpd.clip(l_candidates, df_vias)

    df_l = df_l_cropped

    ## Removendo logradouros pequenos
    df_l = df_l[df_l.length > 6].reset_index()

    df_l = df_l[~(df_l.geometry.type == 'GeometryCollection')]

    pontos = df_l.apply(lambda row: 
                        [row.geometry.interpolate(0, normalized=True),
                        row.geometry.interpolate(1, normalized=True)], axis=1)
    pontos = [item for sublist in pontos for item in sublist]
    df_pontos = gpd.GeoDataFrame(geometry=pontos)

    cruzamento_centroide = [p.centroid for p in df_pontos.buffer(tolerancia_do_cruzamento).unary_union]
    # cruzamento = [nearest_points(p.centroid, df_l.unary_union) for p in df_pontos.buffer(2).unary_union]

    df_cruzamentos = gpd.GeoDataFrame(geometry=cruzamento_centroide)

    df_ct = df_ct.drop(df_ct[(df_ct.geometry.area<200) & (~df_ct.geometry.disjoint(df_cruzamentos.geometry.unary_union))].index, axis=0)

    df_cruzamento_buffer_50 = gpd.GeoDataFrame(geometry=df_cruzamentos.buffer(50))
    df_canteiro_cruzamento = gpd.sjoin(df_ct, df_cruzamento_buffer_50, how="left", op='intersects')

    ps = list(map(lambda x: [x.coords[0][0], x.coords[0][1]], df_cruzamentos.geometry))

    df_voronoi_polygon, df_voronoi_point = voronoi_frames(ps)


    df_quadras_voronoi = gpd.overlay(df_voronoi_polygon, df_vias, how='difference')

    def linhas_curtas(geom, i=0):
        linhas = []
        try:
    #         linhas = []
            for geom in geom:
                linhas.append(LineString(nearest_points(geom, df_cruzamentos.iloc[i].geometry)))
    #         return MultiLineString(linhas)
        except:
            linhas.append(LineString(nearest_points(geom, df_cruzamentos.iloc[i].geometry)))
        return linhas

    linhas_curtas_quadras = [linhas_curtas(row.geometry, i) for i, row in df_quadras_voronoi.iterrows()]
    linhas_curtas_quadras = [item for sublist in linhas_curtas_quadras for item in sublist]
    df_linhas_curtas_quadras = gpd.GeoDataFrame(geometry=linhas_curtas_quadras)

    df_l_out = gpd.overlay(df_l, gpd.GeoDataFrame(geometry=df_cruzamentos.buffer(tolerancia_do_cruzamento)), how='difference')

    # inter_quadra_l = df_linhas_curtas_quadras.apply(lambda row:
    #                                     row.geometry.intersection(df_l_out.unary_union).type == 'Point' or
    #                                     row.geometry.intersection(df_l_out.unary_union).type == 'MultiPoint' , axis=1)

    inter_quadra_l = df_linhas_curtas_quadras.intersects(df_l_out.unary_union)

    linha_cruzamento_canteiro = df_canteiro_cruzamento[~df_canteiro_cruzamento.index_right.isna()].apply(lambda row: 
                                                            LineString(nearest_points(df_cruzamentos.iloc[int(row['index_right'])].geometry, row.geometry)),
                                                            axis=1)

    s = 1.05

    if not linha_cruzamento_canteiro.empty:
        df_linha_cruzamento_canteiro = gpd.GeoDataFrame(geometry=linha_cruzamento_canteiro)

        # inter_canteiro_l = df_linha_cruzamento_canteiro.apply(lambda row: 
        #                                                     row.geometry.intersection(df_l_out.unary_union).type == 'Point' or 
        #                                                     row.geometry.intersection(df_l_out.unary_union).type == 'MultiPoint', 
        #                                                     axis=1)

        inter_canteiro_l = df_linha_cruzamento_canteiro.intersects(df_l_out.unary_union)

        corte_canteiros = df_linha_cruzamento_canteiro[~inter_canteiro_l].apply(lambda x: affinity.scale(x.geometry, xfact=s, yfact=s, origin=x.geometry.coords[0]), axis=1)
        df_corte_canteiros = gpd.GeoDataFrame(geometry=corte_canteiros)
    else:
        df_corte_canteiros = gpd.GeoDataFrame(geometry=[])

    corte_quadras = df_linhas_curtas_quadras[~inter_quadra_l].apply(lambda x: affinity.scale(x.geometry, xfact=s, yfact=s, origin=x.geometry.coords[1]), axis=1)
    df_corte_quadras = gpd.GeoDataFrame(geometry=corte_quadras)

    vias_cutted = unary_union([df_vias.boundary.geometry.unary_union, 
                                            df_ct.boundary.geometry.unary_union])

    linhas_corte = unary_union([df_corte_canteiros[df_corte_canteiros.is_valid].unary_union, 
                               df_corte_quadras[df_corte_quadras.is_valid].unary_union])

    poligonos_de_vias = polygonize(cascaded_union([vias_cutted, linhas_corte, df_s.envelope.boundary.buffer(10).unary_union]))

    df_poligonos_de_vias = gpd.GeoDataFrame(geometry = list(poligonos_de_vias))

    if df_ct.empty:
        df_vias_clipped = df_vias
    else:
        df_vias_clipped = gpd.overlay(df_vias, df_ct, how='difference')

    # df_pvias = df_poligonos_de_vias[df_poligonos_de_vias.geometry.within(df_vias_clipped.buffer(0.1).geometry.unary_union)]
    df_pvias = df_poligonos_de_vias[df_poligonos_de_vias.geometry.intersects(df_vias_clipped.buffer(-0.1).geometry.unary_union)]
    df_pvias.crs = 'EPSG:31983'

    df_pvias.to_file("./resultado/poligono_de_vias.gpkg", layer=f"{list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].lower()}.gpkg", driver="GPKG", OVERWRITE='YES')
    df_pvias.to_file(f"./resultado/{list(df_s.ds_codigo)[0]}_poligono_de_vias_de_{list(df_s.ds_nome)[0].lower().replace(' ', '_')}.gpkg", driver="GPKG")

