import numpy as np
import math
from osgeo import gdal, gdal_array
import os
import matplotlib.pyplot as plt

class CellAreas(object):
    def __init__(self):
        self.b = 6356.7523142
        self.e = 0.0818192
        return

    def import_sample_raster(self, path):
        self.array = gdal_array.LoadFile(path)
        self.nrows = self.array.shape[0]
        self.ncols = self.array.shape[1]
        self.file = gdal.Open(path)
        gtf = list(self.file.GetGeoTransform())
        self.proj = self.file.GetProjection()
        self.gtf = gtf
        self.xcellsize = self.gtf[1]
        self.ycellsize = self.gtf[5]
        self.ul = np.array((gtf[3] , gtf[0]))
        self.cellsize = self.gtf[1]
        self.areas = np.zeros_like(self.array).astype(float)

    def cell_corners_coords(self, cell):
        cell_ul_lat = math.radians(self.ul[0] + self.ycellsize * cell[0])
        cell_ul_lon = math.radians(self.ul[1] + self.xcellsize * cell[1])
        cell_br_lat = cell_ul_lat + math.radians(self.ycellsize)
        cell_br_lon = cell_ul_lon + math.radians(self.xcellsize)
        coords = (cell_ul_lon, cell_br_lon, cell_br_lat, cell_ul_lat)
        return coords

    def w2(self, lat):
        w2 = 1 - (self.e ** 2) * (math.sin(lat) ** 2)
        return w2

    def f(self, lat):
        f = ((self.b ** 2) / 2) * ((math.sin(lat) / (self.w2(lat))) + (1 / (2 * self.e)) * math.log((1 + self.e * math.sin(lat)) / (1 - self.e * math.sin(lat))))
        return f

    def execute(self):
        count = 0
        for i in range(self.nrows):
            j = 1
            count += 1
            cell = np.array((i,j))
            lon1, lon2, lat1, lat2 = self.cell_corners_coords(cell)
            f1 = self.f(lat1)
            f2 = self.f(lat2)
            t = (lon2 - lon1) * (f2 - f1)
            self.areas[i,:] = t

    def save_result(self, path):
        gtf = self.gtf.copy()

        driver = gdal.GetDriverByName('GTiff')

        gtiff = driver.Create(path, self.areas.shape[1], self.areas.shape[0], 1, gdal.GDT_Float64)
        gtiff.SetGeoTransform(gtf)
        gtiff.SetProjection(self.proj)
        gtiff.GetRasterBand(1).WriteArray(self.areas)
        gtiff.FlushCache()
        gtiff = None

test = CellAreas()
test.import_sample_raster('flowdir_rasters/sa_flowdir.tif')
test.execute()
test.save_result('sa_cell_areas.tif')
