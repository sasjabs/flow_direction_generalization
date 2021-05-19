import numpy as np
from osgeo import gdal_array, gdal
import heapq
import os
from geopy.distance import geodesic

class OverallRatio(object):
    def __init__(self):
        self.ratio = 0

    def import_flowdir(self, flowdir_path):
        self.flowdir_path = flowdir_path
        self.flowdir_array = gdal_array.LoadFile(self.flowdir_path)
        file = gdal.Open(self.flowdir_path)
        gtf = list(file.GetGeoTransform())
        proj = file.GetProjection()
        self.nrows = self.flowdir_array.shape[0]
        self.ncols = self.flowdir_array.shape[1]
        self.gtf = gtf
        self.pixelsize = self.gtf[1]
        self.proj = proj

    def calculate_ratio(self):
        orthogonal_count = 0
        diagonal_count = 0
        for i in range(self.nrows):
            for j in range(self.ncols):
                cell = np.array((i,j))
                cell_dir = self.flowdir_array[tuple(cell)]
                if cell_dir in [1, 4, 16, 64]:
                    orthogonal_count += 1
                elif cell_dir in [2, 8, 32, 128]:
                    diagonal_count += 1
        if diagonal_count != 0:
            self.ratio = orthogonal_count / diagonal_count
        else:
            self.ratio = 0

    def execute(self):
        self.calculate_ratio()
        return self.ratio

class WatershedStats(object):
    def __init__(self):
        self.watersheds_mfl = {}
        self.watersheds_ratio = {}
        self.window = {1: np.array((0, 1)),
                       2: np.array((1, 1)),
                       4: np.array((1, 0)),
                       8: np.array((1, -1)),
                       16: np.array((0, -1)),
                       32: np.array((-1, -1)),
                       64: np.array((-1, 0)),
                       128: np.array((-1, 1))}

    def import_rasters(self, flowdir_path, watersheds_path):
        self.flowdir_path = flowdir_path
        self.watersheds_path = watersheds_path
        self.flowdir_array = gdal_array.LoadFile(self.flowdir_path)
        self.watersheds_array = gdal_array.LoadFile(self.watersheds_path)
        file = gdal.Open(self.flowdir_path)
        gtf = list(file.GetGeoTransform())
        proj = file.GetProjection()
        file = None
        self.nrows = self.flowdir_array.shape[0]
        self.ncols = self.flowdir_array.shape[1]
        self.gtf = gtf
        self.xcellsize = self.gtf[1]
        self.ycellsize = self.gtf[5]
        self.ul = np.array((gtf[3] + self.ycellsize/2, gtf[0]+self.xcellsize/2))
        self.proj = proj
        self.watersheds_list = np.unique(self.watersheds_array)
        self.cellsize = self.gtf[1]
        self.flowlen_array = np.zeros_like(self.flowdir_array).astype(float)

    def cell_center_coords(self, cell):
        cell_lat = self.ul[0] + self.ycellsize * cell[0]
        cell_lon = self.ul[1] + self.xcellsize * cell[1]
        coords = np.array((cell_lat, cell_lon))
        return coords

    def watershed_mask(self, watershed_code):
        watershed_mask = (self.watersheds_array == watershed_code).astype(int)
        return watershed_mask

    def watershed_mfl_init(self, watershed_code):
        mask = self.watershed_mask(watershed_code)
        watershed_flowdir = self.flowdir_array * mask
        break_flag = False

        for i in range(self.nrows):
            for j in range(self.ncols):
                cell = np.array((i,j))
                if watershed_flowdir[tuple(cell)] != 0:
                    break_flag = True
                    break
            if break_flag:
                break

        watershed_outlet = self.find_watershed_outlet(cell)

        return watershed_outlet

    def find_watershed_outlet(self, start_cell):
        current_cell = start_cell
        current_dir = self.flowdir_array[tuple(current_cell)]

        while current_dir in [1, 2, 4, 8, 16, 32, 64, 128]:
            next_cell = current_cell + self.window[current_dir]
            current_cell = next_cell
            current_dir = self.flowdir_array[tuple(current_cell)]

        return current_cell

    def calculate_downstream_lengths(self, watershed_outlet):
        q = []
        heapq.heappush(q, watershed_outlet)
        while len(q) > 0:
            center = np.array(heapq.heappop(q))
            for i in self.window.keys():
                current = center + self.window[i]
                if (current[0] < 0) or (current[1] < 0) or (current[0] > self.nrows - 1) or (current[1] > self.ncols - 1):
                    continue
                current_drainage_index = self.flowdir_array[tuple(current)]
                if current_drainage_index not in [1, 2, 4, 8, 16, 32, 64, 128]:
                    continue
                current_rel = self.window[current_drainage_index]
                if np.array_equal(current + current_rel, center):
                    heapq.heappush(q, tuple(current))
                    center_coords = self.cell_center_coords(center)
                    current_coords = self.cell_center_coords(current)
                    dist = geodesic(tuple(current_coords), tuple(center_coords)).kilometers
                    self.flowlen_array[tuple(current)] = self.flowlen_array[tuple(center)] + dist

    def mean_flow_length(self, watershed_code):
        outlet = self.watershed_mfl_init(watershed_code)
        self.calculate_downstream_lengths(outlet)
        watershed_flowlen = self.flowlen_array * ((self.watersheds_array == watershed_code).astype(int))
        if np.count_nonzero(watershed_flowlen) > 0:
            mean_length = np.sum(watershed_flowlen) / np.count_nonzero(watershed_flowlen)
        else: mean_length = 0
        return mean_length

    def calculate_watersheds_mfl(self):
        for watershed_code in self.watersheds_list:
            mfl = self.mean_flow_length(watershed_code)
            self.watersheds_mfl[watershed_code] = mfl

    def watershed_ratio(self, watershed_code):
        mask = self.watershed_mask(watershed_code)
        watershed_flowdir = self.flowdir_array * mask
        orthogonal_sum = 0
        diagonal_sum = 0
        for i in range(self.nrows):
            for j in range(self.ncols):
                current_cell = np.array((i,j))
                cell_dir = watershed_flowdir[tuple(current_cell)]
                if cell_dir in [1, 2, 4, 8, 16, 32, 64, 128]:
                    next_cell = current_cell + self.window[cell_dir]
                    next_cell_coords = self.cell_center_coords(next_cell)
                    current_cell_coords = self.cell_center_coords(current_cell)
                    dist = geodesic(tuple(current_cell_coords), tuple(next_cell_coords)).kilometers
                    if cell_dir in [1, 4, 16, 64]:
                        orthogonal_sum += dist
                    elif cell_dir in [2, 8, 32, 128]:
                        diagonal_sum += dist
        if diagonal_sum > 0:
            ratio = orthogonal_sum / diagonal_sum
        else: ratio = 0
        return ratio

    def export_to_csv(self):
        filename, ext = os.path.splitext(os.path.basename(self.flowdir_path))
        file = open('watershed_stats/'+filename+'_stats.csv', 'w')
        file.write('watershed_code\tmfl_km\tratio\n')
        for watershed_code in self.watersheds_list:
            if watershed_code != 127:
                mfl = self.watersheds_mfl[watershed_code]
                ratio = self.watersheds_ratio[watershed_code]
                line = str(watershed_code)+'\t'+str(mfl)+'\t'+str(ratio)+'\n'
                file.write(line)
        file.close()

    def calculate_watersheds_ratio(self):
        for watershed_code in self.watersheds_list:
            if watershed_code != 127:
                ratio = self.watershed_ratio(watershed_code)
                self.watersheds_ratio[watershed_code] = ratio

    def execute(self):
        self.calculate_watersheds_mfl()
        self.calculate_watersheds_ratio()
        self.export_to_csv()