import pandas as pd
import geopandas as gpd
from osgeo import gdal
from shapely.geometry import Point, Polygon, LineString
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings

warnings.filterwarnings("ignore")

def NTM(in_grid: str, in_lines: str, out_grid: str, out_graph: str, threshold: float):

    # импорт входной сетки
    print('Importing data')
    raster = gdal.Open(in_grid)
    proj = raster.GetProjection()
    gtf = raster.GetGeoTransform()
    arr = raster.GetRasterBand(1).ReadAsArray()
    nodata = raster.GetRasterBand(1).GetNoDataValue()
    grid = arr.copy()

    raster.FlushCache()
    raster = None

    ulx = gtf[0]
    uly = gtf[3]
    dx = gtf[1]
    dy = gtf[5]

    # создание векторной сетки на основе растровой
    cells = []
    irows = []
    icols = []
    for irow in range(grid.shape[0]):
    for icol in range(grid.shape[1]):
        orig_x = ulx + icol * dx
        orig_y = uly + irow * dy
        p0 = Point(orig_x, orig_y)
        p1 = Point(orig_x + dx, orig_y)
        p2 = Point(orig_x + dx, orig_y + dy)
        p3 = Point(orig_x, orig_y + dy)
        cell = Polygon([p0, p1, p2, p3])
        cells.append(cell)
        irows.append(irow)
        icols.append(icol)

    cells = gpd.GeoDataFrame({'row': irows,
                          'col': icols}, 
                         geometry=cells, 
                         crs='epsg:4326')

    # импорт исходных водотоков
    lines = gpd.read_file(in_lines)
    lines['geometry'] = lines['geometry'].translate(1e-4, 1e-4)
    # cells = gpd.read_file('data/Cells.shp').drop(columns = ['Id'])
    cells['cell_id'] = cells.index.copy()

    # пересечение линий с сеткой
    print('Overlaying lines with grid')
    lines_cut = gpd.overlay(lines, cells, how = 'intersection').explode().reset_index(drop=True)
    lines_cut['line_id'] = lines_cut.index.copy()

    # получение начальных и конечных точек линий с нужными атрибутами для связывания графа
    print('Identifying line endpoints')
    startpts = []
    endpts = []
    line_ids = []
    cell_ids = []

    for i, ln in lines_cut.iterrows():
    startpts.append(ln.geometry.boundary[0])
    endpts.append(ln.geometry.boundary[1])
    line_ids.append(ln.line_id)
    cell_ids.append(ln.cell_id)

    startpts = gpd.GeoDataFrame({'lower_line_id': line_ids, 'cell_id': cell_ids}, geometry = startpts, crs = 'epsg:4326')
    endpts = gpd.GeoDataFrame({'upper_line_id': line_ids, 'cell_id': cell_ids}, geometry = endpts, crs = 'epsg:4326')
    startpts['node_id'] = startpts.index.copy()
    endpts['node_id'] = -endpts.index.copy()

    # вместо стандартных операций выборки используется sjoin
    # так как он по умолчанию использует R-tree indexing
    # код выполняется в десятки раз быстрее, чем со стандартными операциями
    print('Identifying border points')
    borderpts = gpd.sjoin(startpts, 
                  gpd.GeoDataFrame(geometry=cells.boundary.buffer(1e-9)), 
                  how='inner').drop_duplicates('node_id').reset_index(drop=True).sort_values('lower_line_id').drop(columns='index_right')

    print('Identifying sources')
    intermed = gpd.sjoin(startpts, 
                 gpd.GeoDataFrame(geometry=endpts.buffer(1e-9)), 
                 how='inner').drop_duplicates('node_id').sort_values('lower_line_id').drop(columns='index_right')

    sourcepts = startpts.loc[list(set(startpts.index).symmetric_difference(set(intermed.index)))].reset_index(drop=True)

    print('Identifying join points')
    intermed1 = gpd.sjoin(startpts, 
                  gpd.GeoDataFrame(geometry=endpts.buffer(1e-9)), 
                  how='inner').drop_duplicates('node_id').sort_values('lower_line_id').drop(columns='index_right')

    intermed2 = gpd.sjoin(intermed1, 
                      gpd.GeoDataFrame(geometry=borderpts.buffer(1e-9)), 
                      how='inner').drop_duplicates('node_id').sort_values('lower_line_id').drop(columns='index_right')

    joinpts = intermed1.loc[list(set(intermed1.index).symmetric_difference(set(intermed2.index)))].reset_index(drop=True)

    print('Identifying mouths')
    intermed3 = gpd.sjoin(endpts, 
                  gpd.GeoDataFrame(geometry=startpts.buffer(1e-9)), 
                  how='inner').drop_duplicates('node_id').sort_values('upper_line_id').drop(columns='index_right')

    mouthpts = endpts.loc[list(set(endpts.index).symmetric_difference(set(intermed3.index)))].reset_index(drop=True).drop(columns=['upper_line_id'])

    borderpts['type'] = 'outlet'
    sourcepts['type'] = 'source'
    joinpts['type'] = 'join'
    mouthpts['type'] = 'mouth'

    # Создание единой таблицы узлов, у каждого узла есть атрибут, связывающий его с вышележащей и нижележащей линией
    print('Linking network nodes')
    endbuff = endpts.copy()
    endbuff.geometry = endbuff.geometry.buffer(1e-6)

    startbuff = startpts.copy()
    startbuff.geometry = startbuff.geometry.buffer(1e-6)

    nodes = pd.concat([borderpts, sourcepts, joinpts, mouthpts], axis=0).reset_index(drop=True).drop(columns='cell_id')

    # здесь для каждой конечной точки при слиянии появляется отдельная запись в таблице
    print('\tJoining...')
    nodes = gpd.sjoin(nodes,
                  endbuff[['geometry', 'upper_line_id', 'cell_id']],
                  how='left').drop(columns=['index_right']).reset_index(drop=True)

    nodes.loc[nodes['type'] == 'source'] = gpd.sjoin(nodes[nodes['type'] == 'source'].drop(columns='cell_id'),
                                                 startbuff[['geometry', 'cell_id']],
                                                 how='left').drop(columns=['index_right']).drop_duplicates('node_id', keep='first')

    # Создание кортежей со ссылками на несколько вышележащих линий для точек слияния
    joins = nodes[nodes['type'] == 'join'].reset_index(drop=True)
    nodes = nodes[nodes['type'] != 'join'].reset_index(drop=True)

    joins_drop = joins.drop_duplicates('node_id').reset_index(drop=True)
    joins_drop.upper_line_id = joins_drop.upper_line_id.astype(object)

    print('\tIdentifying upper lines...')
    for i, row in joins_drop.iterrows():
    nid = row.node_id
    subset = joins[joins.node_id == nid]
    ulids = tuple(subset['upper_line_id'])
    joins_drop.at[i, 'upper_line_id'] = ulids

    nodes = pd.concat([nodes, joins_drop]).reset_index(drop=True)
    nodes['upper_line_id'] = nodes['upper_line_id'].apply(lambda x: (x,) if not isinstance(x, tuple) else x)

    lines_cut['len'] = lines_cut.geometry.length
    nodes = pd.merge(nodes, lines_cut[['line_id', 'len']], 
                 left_on='lower_line_id', right_on='line_id',
                 how='left')

    nodes = nodes.rename(columns={'len': 'down_len'}).drop(columns='line_id')

    # Связывание узлов напрямую, без необходимости каждый раз обращаться к линиям

    print('\tChecking for errors in point types...')
    nodes['up_node'] = None
    nodes['down_node'] = None

    err_indices = []
    for i, row in nodes[nodes['type'] != 'source'].iterrows():
    ulids = row.upper_line_id
    if np.isnan(ulids[0]):
        err_indices.append(i)

    nodes.loc[err_indices, 'type'] = 'source'

    print('\tLinking...')
    cnt = 0
    n_sources = nodes[nodes['type'] != 'source'].shape[0]
    for i, row in nodes[nodes['type'] != 'source'].iterrows():
    if cnt % 1000 == 0:
        print(f'\t\tProcessing point {cnt + 1} out of {n_sources}')
    cnt += 1
    ulids = row.upper_line_id
    up_nodes = []
    for ulid in ulids:
        up_id = nodes[nodes.lower_line_id == ulid].index[0]
        up_nodes.append(up_id)
        nodes.at[up_id, 'down_node'] = i
    nodes.at[i, 'up_node'] = tuple(up_nodes)

    # Расчёт Upstream Length   
    print('Calculating upstream lengths...')
    nodes['up_len'] = 0
    n_sources = nodes[nodes['type'] == 'source'].shape[0]
    cnt = 0
    for i, row in nodes[nodes['type'] == 'source'].iterrows():
    if cnt % 1000 == 0:
        print(f'\tProcessing point {cnt + 1} out of {n_sources}')
    cnt += 1
    ul = 0
    pointer = i
    prev_pointers = []
    while True:
        if nodes.at[pointer, 'type'] == 'join':
            if ul > nodes.at[pointer, 'up_len']:
                nodes.at[pointer, 'up_len'] = ul
                dlen = nodes.at[pointer, 'down_len']
                ul += dlen
                prev_pointers.append(pointer)
                pointer = nodes.at[pointer, 'down_node']
            else:
                break
        elif nodes.at[pointer, 'type'] == 'mouth':
            nodes.at[pointer, 'up_len'] = ul
            break
        else:
            nodes.at[pointer, 'up_len'] = ul
            dlen = nodes.at[pointer, 'down_len']
            ul += dlen
            prev_pointers.append(pointer)
            pointer = nodes.at[pointer, 'down_node']

    # Определение выходных точек ячеек            
    print('Identifying cell exit points')
    cells['exit_node'] = None
    notna_cells = set(nodes['cell_id'])
    for i, row in cells.iterrows():
    if i not in notna_cells:
        continue

    exit_node = nodes[nodes['cell_id'] == i].sort_values('up_len', ascending=False).index[0]
    cells.at[i, 'exit_node'] = exit_node

    # Определение направления стока из ячейки    
    print('Calculating cells flow direction')

    cells['down_cell'] = None
    cells['fdir'] = None

    window = {0: (0, 0),
          1: (0, 1),
          2: (1, 1),
          4: (1, 0),
          8: (1, -1),
          16: (0, -1),
          32: (-1, -1),
          64: (-1, 0),
          128: (-1, 1)}

    rev_window = {(0, 0): 0,
              (0, 1): 1,
              (1, 1): 2,
              (1, 0): 4,
              (1, -1): 8,
              (0, -1): 16,
              (-1, -1): 32,
              (-1, 0): 64,
              (-1, 1): 128}

    n_cells = cells[cells['exit_node'].notna()].shape[0]
    cnt = 0
    for i, row in cells[cells['exit_node'].notna()].iterrows():

    if cnt % 100 == 0:
        print(f'\tProcessing cell {cnt + 1} out of {n_cells}')
    cnt += 1

    pointer = row.exit_node
    start_ul = nodes.at[pointer, 'up_len']

    while True:
        if nodes.at[pointer, 'type'] == 'mouth':
            break

        pointer = nodes.at[pointer, 'down_node']

        if nodes.at[pointer, 'type'] == 'outlet':
            cur_ul = nodes.at[pointer, 'up_len']
            stop = 0

            if (cur_ul - start_ul) > threshold:
                stop = 1
                break
            elif stop == 1:
                break

    down_cell = nodes.at[pointer, 'cell_id']
    down_cell_row = cells.at[down_cell, 'row']
    down_cell_col = cells.at[down_cell, 'col']

    curr_cell_row = cells.at[i, 'row']
    curr_cell_col = cells.at[i, 'col']

    dx = down_cell_col - curr_cell_col
    dy = down_cell_row - curr_cell_row

    if abs(dx) > 1 or abs(dy) > 1:
        pointer = row.exit_node

        while True:
            if nodes.at[pointer, 'type'] == 'mouth':
                break

            pointer = nodes.at[pointer, 'down_node']

            if nodes.at[pointer, 'type'] == 'outlet':
                break

        down_cell = nodes.at[pointer, 'cell_id']
        down_cell_row = cells.at[down_cell, 'row']
        down_cell_col = cells.at[down_cell, 'col']

        curr_cell_row = cells.at[i, 'row']
        curr_cell_col = cells.at[i, 'col']

        dx = down_cell_col - curr_cell_col
        dy = down_cell_row - curr_cell_row

    fdir = rev_window[(dy, dx)]
    grid[curr_cell_row, curr_cell_col] = fdir

    cells.at[i, 'fdir'] = fdir

    if down_cell == i:
        cells.at[i, 'down_cell'] = 0
    else:
        cells.at[i, 'down_cell'] = down_cell

    # Выявление ячеек, которые нужно перетрассировать из-за наличия пересечений
    print('Retracing cells')
    retrace_cells = []
    for i in range(grid.shape[0] - 1):
    for j in range(grid.shape[1] - 1):
        # print(grid[i, j], grid[i, j+1])
        # print(grid[i+1, j], grid[i+1, j+1])
        # print('\n')
        if grid[i + 1, j] == 128 and grid[i + 1, j + 1] == 32:

            cell1 = cells[(cells['row'] == i + 1) & 
                          (cells['col'] == j)].index[0]

            cell2 = cells[(cells['row'] == i + 1) & 
                          (cells['col'] == j + 1)].index[0]

            retrace_cells.extend((cell1, cell2))

        if grid[i, j] == 2 and grid[i, j + 1] == 8:

            cell1 = cells[(cells['row'] == i) & 
                          (cells['col'] == j)].index[0]

            cell2 = cells[(cells['row'] == i) & 
                          (cells['col'] == j + 1)].index[0]

            retrace_cells.extend((cell1, cell2))

        if grid[i, j + 1] == 8 and grid[i + 1, j + 1] == 32:

            cell1 = cells[(cells['row'] == i) & 
                          (cells['col'] == j + 1)].index[0]

            cell2 = cells[(cells['row'] == i + 1) & 
                          (cells['col'] == j + 1)].index[0]

            retrace_cells.extend((cell1, cell2))

        if grid[i, j] == 2 and grid[i + 1, j] == 128:

            cell1 = cells[(cells['row'] == i) & 
                          (cells['col'] == j)].index[0]

            cell2 = cells[(cells['row'] == i + 1) & 
                          (cells['col'] == j)].index[0]

            retrace_cells.extend((cell1, cell2))

    # Исправление пересечений    
    print('Fixing intersections')
    for icell in retrace_cells:
    pointer = cells.at[icell, 'exit_node']
    while True:
        if nodes.at[pointer, 'type'] == 'mouth':
            break

        pointer = nodes.at[pointer, 'down_node']

        if nodes.at[pointer, 'type'] == 'outlet':
            break

    down_cell = nodes.at[pointer, 'cell_id']
    down_cell_row = cells.at[down_cell, 'row']
    down_cell_col = cells.at[down_cell, 'col']

    curr_cell_row = cells.at[icell, 'row']
    curr_cell_col = cells.at[icell, 'col']

    dx = down_cell_col - curr_cell_col
    dy = down_cell_row - curr_cell_row

    fdir = rev_window[(dy, dx)]
    # print(dy, dx, fdir)
    grid[curr_cell_row, curr_cell_col] = fdir

    cells.at[icell, 'fdir'] = fdir

    if down_cell == i:
        cells.at[icell, 'down_cell'] = 0
    else:
        cells.at[icell, 'down_cell'] = down_cell

    print('Creating vector graph')
    graph = cells.copy()
    graph.geometry = graph.geometry.centroid
    lines = []
    for i, row in graph[graph['down_cell'].notna()].iterrows():
    if row.down_cell == 0:
        continue
    start = row.geometry
    end = graph.at[row.down_cell, 'geometry']
    lines.append(LineString([start, end]))

    # Сохранение результата    
    print('Saving the result')
    lines_graph = gpd.GeoDataFrame(geometry=lines, crs='epsg:4326')

    driver = gdal.GetDriverByName('GTiff')
    gtiff = driver.Create(out_grid, grid.shape[1], grid.shape[0], 1, gdal.GDT_Int16)
    gtiff.SetGeoTransform(gtf)
    gtiff.SetProjection(proj)
    gtiff.GetRasterBand(1).SetNoDataValue(nodata)
    gtiff.GetRasterBand(1).WriteArray(grid)
    gtiff.FlushCache()
    gtiff = None

    lines_graph.to_file(out_graph, index=False)