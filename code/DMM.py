import numpy as np
from osgeo import gdal_array, gdal  

class DMM(object):
    def __init__(self, k: int):
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
        self.cells = np.full(
            (self.flowacc_array.shape[0] // self.k - 1, self.flowacc_array.shape[1] // self.k - 1), 255, dtype=np.uint8)
        self.mrows = self.cells.shape[0]
        self.mcols = self.cells.shape[1]
        self.null_mask = (self.flowacc_array != self.ndv).astype(int)

    def save_result(self, path):
        gtf = self.gtf.copy()
        k = self.k
        gtf[1] = gtf[1] * k
        gtf[5] = gtf[5] * k

        driver = gdal.GetDriverByName('GTiff')

        gtiff = driver.Create(path, self.cells.shape[1], self.cells.shape[0], 1, gdal.GDT_Int32)
        gtiff.SetGeoTransform(gtf)
        gtiff.SetProjection(self.proj)
        gtiff.GetRasterBand(1).SetNoDataValue(255)
        gtiff.GetRasterBand(1).WriteArray(self.cells)
        gtiff.FlushCache()
        gtiff = None

    def stack_rows_and_cols(self, array):
        nrows = array.shape[0]
        ncols = array.shape[1]
        rows_rem = nrows % self.k
        cols_rem = ncols % self.k

        array_add = np.concatenate((array, np.zeros([nrows, cols_rem], dtype=int)), axis=1)
        array_add = np.concatenate((array_add, np.zeros([rows_rem, ncols + cols_rem], dtype=int)), axis=0)

        array_add = np.pad(array_add, pad_width=self.k // 2, mode='constant', constant_values=0)

        return array_add

    def find_max(self, pixels):  # находит индексы максимального или минимального элемента массива
        result = np.where(pixels == np.amax(pixels))
        coords = np.squeeze(np.array(list(zip(result[0], result[1]))))
        return coords

    def find_largest_drainage_area_pixel(self,
                                         cell,
                                         displaced=False):  # находит внутри ячейки пиксель с наибольшей водосборной площадью
        subset = self.get_pixels_in_cell(self.flowacc_array, cell, displaced)
        pixel = self.find_max(subset)
        if len(pixel.shape) > 1:
            pixel = pixel[0]
        pixel_global = self.cell_coords_to_global(pixel, cell, displaced)
        return pixel_global

    def get_pixels_in_cell(self, array, cell, displaced=False):  # возвращает пиксели в пределах ячейки
        i, j = cell
        k = self.k
        if displaced:
            d = 0
        else:
            d = k // 2
        left = i * k + d
        right = i * k + k + d
        top = j * k + d
        bottom = j * k + k + d
        pixels = array[left:right, top:bottom].copy()
        return pixels

    def identify_cell_by_pixel(self, pixel, displaced=False):
        if displaced:
            d = 0
        else:
            d = self.k // 2
        cell = np.array(((pixel[0] - d) // self.k, (pixel[1] - d) // self.k))
        return cell

    def cell_coords_to_global(self, pixel_coords, cell,
                              displaced=False):  # пересчитывает координаты пикселя в пределах ячейки в "глобальные"
        x, y = cell
        m, n = pixel_coords
        k = self.k
        if displaced:
            d = 0
        else:
            d = self.k // 2
        global_coords = np.array([x * k + m + d, y * k + n + d])
        return global_coords

    def define_direction(self, diff):
        dir = 0

        if np.array_equal(diff, np.array((0, 1))):
            dir = 1
        elif np.array_equal(diff, np.array((1, 1))):
            dir = 2
        elif np.array_equal(diff, np.array((1, 0))):
            dir = 4
        elif np.array_equal(diff, np.array((1, -1))):
            dir = 8
        elif np.array_equal(diff, np.array((0, -1))):
            dir = 16
        elif np.array_equal(diff, np.array((-1, -1))):
            dir = 32
        elif np.array_equal(diff, np.array((-1, 0))):
            dir = 64
        elif np.array_equal(diff, np.array((-1, 1))):
            dir = 128

        return dir

    def check_null_cell(self, cell):
        pixels = self.get_pixels_in_cell(self.null_mask, cell)
        if np.max(pixels) == 0:
            return True
        else:
            return False

    def check_cell_out_of_bounds(self, cell):
        if np.min(cell) < 0 or cell[0] > self.mrows - 1 or cell[1] > self.mcols - 1:
            return True
        else:
            return False

    def assign_cell_directions(self):
        for i in list(range(self.mrows)):
            for j in list(range(self.mcols)):
                cell = np.array((i, j))
                discharge_cell = cell.copy()
                if not self.check_null_cell(cell):
                    cell_outlet = self.find_largest_drainage_area_pixel(cell, displaced=False)
                    displaced_cell = self.identify_cell_by_pixel(cell_outlet, displaced=True)
                    displaced_outlet = self.find_largest_drainage_area_pixel(displaced_cell, displaced=True)
                    discharge_cell = self.identify_cell_by_pixel(displaced_outlet, displaced=False)
                    if np.array_equal(cell, discharge_cell):
                        max_flowacc = 0
                        max_flowacc_cell = cell
                        for l in self.window.keys():
                            cell_in_window = cell + self.window[l]
                            if not self.check_cell_out_of_bounds(cell_in_window) and not self.check_null_cell(cell_in_window):
                                pixels = self.get_pixels_in_cell(self.flowacc_array, cell_in_window)
                                cell_flowacc = np.max(pixels)
                                if cell_flowacc > max_flowacc:
                                    max_flowacc = cell_flowacc
                                    max_flowacc_cell = cell_in_window
                        discharge_cell = max_flowacc_cell

                diff = discharge_cell - cell
                cell_dir = self.define_direction(diff)
                self.cells[tuple(cell)] = cell_dir

    def fix_counter_flows(self):
        for i in range(self.mrows - 1):
            for j in range(self.mcols - 1):
                pixels_c = self.get_pixels_in_cell(self.flowacc_array, np.array((i,j)), displaced=False)
                pixels_r = self.get_pixels_in_cell(self.flowacc_array, np.array((i,j+1)), displaced=False)
                pixels_b = self.get_pixels_in_cell(self.flowacc_array, np.array((i+1,j)), displaced=False)
                pixels_br = self.get_pixels_in_cell(self.flowacc_array, np.array((i+1,j+1)), displaced=False)
                pixels_tr = self.get_pixels_in_cell(self.flowacc_array, np.array((i-1,j+1)), displaced=False)
                if self.cells[i, j] == 1 and self.cells[i, j+1] == 16:
                    if np.max(pixels_c) > np.max(pixels_r):
                        self.cells[i, j] = 0
                    else:
                        self.cells[i, j+1] = 0
                if self.cells[i, j] == 4 and self.cells[i+1, j] == 64:
                    if np.max(pixels_c) > np.max(pixels_b):
                        self.cells[i, j] = 0
                    else:
                        self.cells[i + 1, j] = 0
                if self.cells[i, j] == 2 and self.cells[i+1, j+1] == 32:
                    if np.max(pixels_c) > np.max(pixels_br):
                        self.cells[i, j] = 0
                    else:
                        self.cells[i+1, j+1] = 0
                if i >= 1 and self.cells[i, j] == 128 and self.cells[i-1, j+1] == 8:
                    if np.max(pixels_c) > np.max(pixels_tr):
                        self.cells[i, j] = 0
                    else:
                        self.cells[i-1, j+1] = 0

    def fix_intersections(self):
        for i in range(self.mrows - 1):
            for j in range(self.mcols - 1):
                pixels_c = self.get_pixels_in_cell(self.flowacc_array, np.array((i, j)), displaced=False)
                pixels_r = self.get_pixels_in_cell(self.flowacc_array, np.array((i, j + 1)), displaced=False)
                pixels_b = self.get_pixels_in_cell(self.flowacc_array, np.array((i + 1, j)), displaced=False)
                pixels_br = self.get_pixels_in_cell(self.flowacc_array, np.array((i + 1, j + 1)), displaced=False)
                if self.cells[i + 1, j] == 128 and self.cells[i + 1, j + 1] == 32:
                    if np.max(pixels_c) > np.max(pixels_r):
                        self.cells[i + 1, j] = 64
                    else:
                        self.cells[i + 1, j + 1] = 64

                if self.cells[i, j] == 2 and self.cells[i, j + 1] == 8:
                    if np.max(pixels_br) > np.max(pixels_b):
                        self.cells[i, j + 1] = 4
                    else:
                        self.cells[i, j] = 4

                if self.cells[i, j + 1] == 8 and self.cells[i + 1, j + 1] == 32:
                    if np.max(pixels_c) > np.max(pixels_b):
                        self.cells[i, j + 1] = 16
                    else:
                        self.cells[i + 1, j + 1] = 16

                if self.cells[i, j] == 2 and self.cells[i + 1, j] == 128:
                    if np.max(pixels_r) > np.max(pixels_br):
                        self.cells[i, j] = 1
                    else:
                        self.cells[i + 1, j] = 1

    def execute(self):
        self.assign_cell_directions()
        self.fix_counter_flows()
        self.fix_intersections()