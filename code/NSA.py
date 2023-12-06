import numpy as np
from osgeo import gdal_array, gdal  

class NSA(object):
    def __init__(self, k):
        self.k = k

        self.window = {1: np.array((0, 1)),
                       2: np.array((1, 1)),
                       4: np.array((1, 0)),
                       8: np.array((1, -1)),
                       16: np.array((0, -1)),
                       32: np.array((-1, -1)),
                       64: np.array((-1, 0)),
                       128: np.array((-1, 1))}

    def import_raster(self, path):  # функция импорта файла, возвращает массив, GeoTransform и Projection.
        array = gdal_array.LoadFile(path)
        file = gdal.Open(path)
        gtf = list(file.GetGeoTransform())
        proj = file.GetProjection()
        ndv = file.GetRasterBand(1).GetNoDataValue()
        return array, gtf, proj, ndv

    def import_flowacc(self, flowacc_path):
        flowacc_array, gtf, proj, ndv = self.import_raster(flowacc_path)
        self.flowacc_array = self.stack_rows_and_cols(flowacc_array)
        self.nrows = self.flowacc_array.shape[0]
        self.ncols = self.flowacc_array.shape[1]
        self.gtf = gtf
        self.pixelsize = self.gtf[1]
        self.proj = proj
        self.ndv = ndv
        self.cells = np.full((self.flowacc_array.shape[0] // self.k, self.flowacc_array.shape[1] // self.k), 255, dtype=np.uint8)
        self.flowacc_aggr = self.cells.copy()
        self.krows = self.cells.shape[0]
        self.kcols = self.cells.shape[1]

    def stack_rows_and_cols(self, array):
        nrows = array.shape[0]
        ncols = array.shape[1]
        rows_rem = nrows % self.k
        cols_rem = ncols % self.k

        array_add = np.concatenate((array, np.zeros([nrows, cols_rem], dtype=int)), axis=1)
        array_add = np.concatenate((array_add, np.zeros([rows_rem, ncols+cols_rem], dtype=int)), axis=0)

        return array_add

    def aggregate(self): # функция агрегации массива, k - размер ячейки агрегированного массива в ячейках исходного
        k = self.k

        for i in range(self.krows):
            for j in range(self.kcols):
                pixels = self.flowacc_array[i*k:i*k+k, j*k:j*k+k]
                self.flowacc_aggr[i,j] = np.max(pixels)

    def check_cell_out_of_bounds(self, cell):
        if np.min(cell) < 0 or cell[0] > self.krows - 1 or cell[1] > self.kcols - 1:
            return True
        else:
            return False

    def flowdir(self): # расчет направлений стока при помощи очереди с приоритетом
        for i in range(self.krows):
            for j in range(self.kcols):
                cell = np.array((i,j))
                cell_flowacc = self.flowacc_aggr[tuple(cell)]
                if cell_flowacc > 0:
                    max_n_grad = 0
                    max_grad_dir = 0
                    for dir in self.window.keys():
                        n_cell = cell + self.window[dir]
                        if not self.check_cell_out_of_bounds(n_cell):
                            n_cell_flowacc = self.flowacc_aggr[tuple(n_cell)]
                            if dir in [1, 4, 16, 64]:
                                len = 1
                            else:
                                len = 2**(0.5)
                            n_cell_grad = (n_cell_flowacc - cell_flowacc) / len
                            if n_cell_grad > max_n_grad:
                                max_n_grad = n_cell_grad
                                max_grad_dir = dir
                    self.cells[tuple(cell)] = max_grad_dir


    def save_result(self, path):
        gtf = self.gtf.copy()
        k = self.k
        gtf[1] = gtf[1] * k
        gtf[5] = gtf[5] * k

        driver = gdal.GetDriverByName('GTiff')

        gtiff = driver.Create(path, self.cells.shape[1], self.cells.shape[0], 1, gdal.GDT_Int32)
        gtiff.SetGeoTransform(gtf)
        gtiff.SetProjection(self.proj)
        gtiff.GetRasterBand(1).WriteArray(self.cells)
        gtiff.GetRasterBand(1).SetNoDataValue(255)
        gtiff.FlushCache()
        gtiff = None

    def execute(self):
        self.aggregate()
        self.flowdir()
