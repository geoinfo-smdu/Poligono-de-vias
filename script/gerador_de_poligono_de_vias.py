#!/usr/bin/env geo_env
# coding: utf-8

# # Gerando polígonos de vias por setor fiscal município de São Paulo
# 
# Esse notebook tem a intenção de ir um pouco além do último Notebook no sentido de avançar para um processamento viável para a publicação de uma nova camada no [Geosampa](http://geosampa.prefeitura.sp.gov.br/PaginasPublicas/_SBC.aspx) 
# 
# Para isso vamos utilizar as seguintes camadas, recortadas por um setor fiscal definido para testes
# 
# * [Logradouros](http://geosampa.prefeitura.sp.gov.br/PaginasPublicas/downloadIfr.aspx?orig=DownloadCamadas&arq=05_Sistema%20Vi%E1rio%5C%5CLogradouro%5C%5CShapefile%5C%5CSIRGAS_SHP_logradouronbl&arqTipo=Shapefile)
# * [Quadras Fiscais](http://geosampa.prefeitura.sp.gov.br/PaginasPublicas/downloadIfr.aspx?orig=DownloadCamadas&arq=11_Cadastro%5C%5CQuadra%5C%5CShapefile%5C%5CSIRGAS_SHP_quadraMDSF&arqTipo=Shapefile)
# * [Quadras Viárias](http://geosampa.prefeitura.sp.gov.br/PaginasPublicas/downloadIfr.aspx?orig=DownloadCamadas&arq=05_Sistema%20Vi%E1rio%5C%5CQuadra%20Viaria%5C%5CShapefile%5C%5CSIRGAS_SHP_quadraviariaed&arqTipo=Shapefile)
# 
# Para nossos primeiros experimentos vamos utilizar o setor fiscal de número 037:
# 
# * [Setores Fiscais](http://geosampa.prefeitura.sp.gov.br/PaginasPublicas/downloadIfr.aspx?orig=DownloadCamadas&arq=11_Cadastro%5C%5CSetor%5C%5CShapefile%5C%5CSIRGAS_SHP_setorfiscal&arqTipo=Shapefile)
# 
# Os arquivos baixados nos links acima devem ser descompactados e colocados na pas `gis` deste repositório para funcioarem de acordo com os scripts a seguir.
# 

# ## Trabalhando com o Geopandas
# 
# O Geopandas é um projeto de código aberto escrito em Python, possui bastante maturidade, eficácia e performance para trabalhar com grandes quantidade de dados. Ele usa o Pandas, NumPy e Shapely para trabalhar com dados georeferenciados e portanto achamos oportuno usa-lo.

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from datetime import datetime

warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)

plt.rcParams['figure.figsize'] = (40, 30)

# df_qf_s = gpd.read_file(f'gis/SIRGAS_SHP_quadraMDSF.shp')
df_qv_s = gpd.read_file(f'gis/SIRGAS_SHP_quadraviariaed_polygon.shp')
df_l_s = gpd.read_file(f'gis/SIRGAS_SHP_logradouronbl.shp')
df_s_s = gpd.read_file(f'gis/SIRGAS_SHP_distrito_polygon.shp')

for i in df_s_s.itertuples():
    df_s = df_s_s[df_s_s.ds_codigo == i.ds_codigo]

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print("****** início: ", dt_string)	
    print(f"Executando {list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].lower()}")
    #df_s = df_s[df_s.ds_codigo == '1']
    # df_s = df_s.query("st_codigo == '037'")

    df_s_sindex = df_s.sindex
    df_l_sindex = df_l_s.sindex

    # Get the bounding box coordinates of the Polygon as a list
    bounds = list(df_s.bounds.values[0])

    # Get the indices of the Points that are likely to be inside the bounding box of the given Polygon
    l_candidate_idx = list(df_l_sindex.intersection(bounds))
    l_candidates = df_l_s.loc[l_candidate_idx]

    df_qv_sindex = df_qv_s.sindex

    # Get the bounding box coordinates of the Polygon as a list
    bounds = list(df_s.bounds.values[0])

    # Get the indices of the Points that are likely to be inside the bounding box of the given Polygon
    qv_candidate_idx = list(df_qv_sindex.intersection(bounds))
    qv_candidates = df_qv_s.loc[qv_candidate_idx]

    df_vias = gpd.overlay(df_s, qv_candidates, how='difference')

    df_l_cropped = gpd.clip(l_candidates, df_vias)

    df_l = df_l_cropped

    df_vias_buffered = df_vias.buffer(0.1)

    # Obtendo os pontos finais e iniciais de cada logradouro podemos saber onde as ruas começam e terminam
    gdf_cruzamentos_inicio = df_l[df_l['geometry'].type == 'LineString']['geometry'].map(lambda x: (x.interpolate(0, normalized=True)))
    gdf_cruzamentos_final = df_l[df_l['geometry'].type == 'LineString']['geometry'].map(lambda x: (x.interpolate(1, normalized=True)))
    gdf_cruzamentos =  pd.concat([gdf_cruzamentos_inicio, gdf_cruzamentos_final])
    gdf_c = gdf_cruzamentos


    # Como alguns pontos acabam ficando muito próximos, seria interessante agregar todos os pontos do cruzamento.
    cruzamentos = [p.centroid for p in gdf_c.buffer(4).unary_union]
    df_cruzamentos = gpd.GeoDataFrame(geometry=cruzamentos)

    df_pts = gpd.GeoDataFrame(geometry=df_cruzamentos.buffer(10).boundary.intersection(df_l.unary_union))
    # filtrando somente os multipontos
    df_cruzamentos = df_cruzamentos[df_pts.geometry.type == 'MultiPoint']
    df_pts = df_pts[df_pts.geometry.type == 'MultiPoint']

    df_cruzamentos['id'] = df_cruzamentos.index.values
    df_pts['id'] = df_pts.index.values
    df_pt_pts = df_cruzamentos.merge(df_pts, on='id', suffixes=('_cruzamento', '_pt'))

    import numpy as np

    def calc_angulos(i):
        c_x, c_y = i.geometry_cruzamento.coords[0]
        coords = list(map(lambda x: x.coords[0], i.geometry_pt.geoms))
        p_x = np.array(list(map(lambda x: x[0], coords)))
        p_y = np.array(list(map(lambda x: x[1], coords)))
        delta_x, delta_y = [p_x - c_x, p_y - c_y]
        # Colocar os pontos em ordem de angulo
        angle = np.sort(np.arctan2(delta_x, delta_y) - 0.5 * np.pi)
        angle_normalized = angle / (np.pi * 2) + 0.75
        angle_plus = np.append(angle_normalized, 1 + angle_normalized[0:1])
        points = angle_normalized + (np.diff(angle_plus) / 2)
        points[-1] = points[-1] - 1
        return list(points)

    angles = list(map(lambda x: calc_angulos(x), df_pt_pts.itertuples()))

    df_bounds = gpd.GeoDataFrame(geometry=df_cruzamentos.buffer(10).boundary)
    df_bounds['angles'] = angles

    from shapely.geometry import MultiPoint

    inter_points = list(map(
        lambda x: MultiPoint(list(map(
            lambda y: x.geometry.interpolate(y - 0.75, normalized=True), 
            x.angles))), 
        df_bounds.itertuples()))

    df_inter_points = gpd.GeoDataFrame(geometry=inter_points)

    # ## Criando linhas para 'cortar' os polígonos
    # 
    # Agora que temos que cortar o polígono precisamos de linhas para realizar esse fatiamento. 
    # Primeiramente vamos criar uma intersecção dos buffers dos nós com o traçado dos logradouros

    from shapely.geometry import MultiLineString, LineString
    from shapely.affinity import scale

    cut_lines = list(map(lambda x:
            scale(
                MultiLineString(list(map(lambda y: 
                    LineString([df_cruzamentos.geometry.iloc[x.Index], y]), 
                    x.geometry.geoms))),
                xfact = 2.0,
                yfact = 2.0,
                origin = df_cruzamentos.geometry.iloc[x.Index])
            , df_inter_points.itertuples()))

    df_cut_lines = gpd.GeoDataFrame(geometry=cut_lines)

    df_cut_lines = gpd.clip(df_cut_lines, df_vias_buffered)

    lp = list(map(lambda x: MultiPoint([x.geometry.interpolate(0.5, normalized=True),
                                x.geometry.interpolate(0.001, normalized=True),
                                x.geometry.interpolate(0.999, normalized=True)]),
            df_l[df_l['geometry'].type == 'LineString'].itertuples()))

    df_lp = gpd.GeoDataFrame(geometry=lp)

    from shapely.ops import cascaded_union, polygonize
    df_vias_cutted = polygonize(cascaded_union([df_cut_lines.geometry.unary_union, list(df_vias.geometry.boundary)[0]]))
    df_poligonos_de_vias = gpd.GeoDataFrame(geometry = list(df_vias_cutted))
    df_poligonos_de_vias.sindex
    # df_pvias = df_poligonos_de_vias[df_poligonos_de_vias.geometry.crosses(df_lp.geometry.unary_union)]
    df_pvias = df_poligonos_de_vias[df_poligonos_de_vias.geometry.within(df_vias_buffered.geometry.unary_union)]

    # ## Salvando o resultado
    # 
    # Agora que temos o polígono de cada via, podemos salvar o resultado para seguir com alguma outra análise como:
    # - Calcular a largura mínima máxima e média de cada via

    df_pvias.to_file("./resultado/poligono_de_vias.gpkg", layer=f"{list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].capitalize()}", driver="GPKG")
    df_pvias.to_file(f"./resultado/poligono_de_vias_{list(df_s.ds_codigo)[0]}_{'_'.join(list(df_s.ds_nome)[0].lower().split())}.gpkg", layer=f"{list(df_s.ds_codigo)[0]} - poligono de vias de {list(df_s.ds_nome)[0].capitalize()}", driver="GPKG")