import numpy as np
from osgeo import gdal_array, gdal  


class COTAT(object):
    def __init__(self, k: int, area_threshold: float):
        self.k = k
        self.area_threshold = area_threshold

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

    def import_flowdir(self, flowdir_path):
        flowdir_array, gtf, proj, ndv = self.import_raster(flowdir_path)
        # print(flowdir_array.dtype)
        # self.flowdir_array = flowdir_array
        self.flowdir_array = self.stack_rows_and_cols(flowdir_array)
        del flowdir_array
        self.nrows = self.flowdir_array.shape[0]
        self.ncols = self.flowdir_array.shape[1]
        self.gtf = gtf
        self.pixelsize = self.gtf[1]
        self.proj = proj
        self.flowdir_ndv = ndv
        self.cells = np.full((self.flowdir_array.shape[0] // self.k, self.flowdir_array.shape[1] // self.k), 255, dtype=np.uint8)
        self.mrows = self.cells.shape[0]
        self.mcols = self.cells.shape[1]

    def create_null_mask(self):
        self.null_mask = (self.flowacc_array != self.flowacc_ndv).astype(int)

    def import_flowacc(self, flowacc_path):
        flowacc_array, gtf, proj, ndv = self.import_raster(flowacc_path)
        # print(flowacc_array.dtype)
        # self.flowacc_array = flowacc_array
        self.flowacc_array = self.stack_rows_and_cols(flowacc_array)
        self.flowacc_ndv = ndv
        self.flowacc_array = np.where(self.flowacc_array == ndv, 0, self.flowacc_array)
        del flowacc_array

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
        
        if rows_rem == 0 and cols_rem == 0:
            
            return array
        
        else:
            
            return np.concatenate((np.concatenate((array, np.zeros([nrows, cols_rem], dtype=np.uint32)), axis=1), np.zeros([rows_rem, ncols + cols_rem], dtype=np.uint32)), axis=0)

    def find_max(self, pixels):  # находит индексы максимального или минимального элемента массива
        result = np.where(pixels == np.amax(pixels))
        coords = np.squeeze(np.array(list(zip(result[0], result[1]))))
        return coords

    def find_largest_drainage_area_pixel(self, cell):  # находит внутри ячейки пиксель с наибольшей водосборной площадью
        subset = self.get_pixels_in_cell(self.flowacc_array, cell)
        pixel = self.find_max(subset)
        if len(pixel.shape) > 1:
            pixel = pixel[0]
        pixel_global = self.cell_coords_to_global(pixel, cell)
        return pixel_global

    def get_pixels_in_cell(self, array, cell):  # возвращает пиксели в пределах ячейки
        i, j = cell
        k = self.k
        left = i * k
        right = i * k + k
        top = j * k
        bottom = j * k + k
        pixels = array[left:right, top:bottom].copy()
        return pixels

    def identify_cell_by_pixel(self, pixel):
        cell = np.array((pixel[0] // self.k, pixel[1] // self.k))
        return cell

    def find_cell_outlet_by_mask(self, cell):
        pixels = self.get_pixels_in_cell(self.cell_outlets, cell)
        result = np.where(pixels == 1)
        outlet = np.squeeze(np.array(list(zip(result[0], result[1]))))
        outlet_global = self.cell_coords_to_global(outlet, cell)
        return outlet_global

    def cell_coords_to_global(self, pixel_coords,
                              cell):  # пересчитывает координаты пикселя в пределах ячейки в "глобальные"
        x, y = cell
        m, n = pixel_coords
        k = self.k
        global_coords = np.array([x * k + m, y * k + n])
        return global_coords

    def border_mask(self):  # создаёт маску для обхода краевых точек
        k = self.k
        mask = np.zeros((k, k), dtype=int)
        mask[0, :] = 1
        mask[:, 0] = 1
        mask[k - 1, :] = 1
        mask[:, k - 1] = 1
        return mask

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

    def check_same_pixel(self, pixel1, pixel2):
        if pixel1[0] == pixel2[0] and pixel1[1] == pixel2[1]:
            return True
        else:
            return False

    def check_cell_change(self, prev_pixel, current_pixel):
        prev_cell = self.identify_cell_by_pixel(prev_pixel)
        curr_cell = self.identify_cell_by_pixel(current_pixel)
        cell_change = not np.array_equal(prev_cell, curr_cell)
        return cell_change

    def check_null_pixel(self, pixel):
        if self.null_mask[tuple(pixel)] == 0:
            return True
        else:
            return False

    def check_null_cell(self, cell):
        pixels = self.get_pixels_in_cell(self.null_mask, cell)
        if np.max(pixels) == 0:
            return True
        else:
            return False

    def check_pixel_out_of_bounds(self, pixel):
        if np.min(pixel) < 0 or pixel[0] > self.nrows - 1 or pixel[1] > self.ncols - 1 or self.check_null_pixel(pixel):
            return True
        else:
            return False

    def check_cell_out_of_bounds(self, center_cell, neighbor_cell):
        if np.max(np.abs(neighbor_cell - center_cell)) > 1:
            return True
        else:
            return False

    def cell_outlet_pixel(self, cell):
        outlet_pixel = self.find_largest_drainage_area_pixel(cell)
        return outlet_pixel

    def assign_outlet_pixels(self):
        self.cell_outlets = np.zeros_like(self.flowdir_array)
        ncells = self.cells.shape[0] * self.cells.shape[1]
        cnt = 0
        interval = ncells // 100
        for i in range(self.cells.shape[0]):
            for j in range(self.cells.shape[1]):
                if cnt % interval == 0:
                    print(f'\t\tProcessing cell {cnt + 1} out of {ncells}')
                cnt += 1
                cell = np.array((i, j))
                if not self.check_null_cell(cell):
                    outlet = self.cell_outlet_pixel(cell)
                    self.cell_outlets[tuple(outlet)] = 1

    def assign_cell_direction(self, cell):
        window = self.window.copy()
        current_pixel = self.find_cell_outlet_by_mask(cell)
        latest_outlet = current_pixel.copy()
        current_area = self.flowacc_array[tuple(current_pixel)]
        init_area = current_area
        area_diff = 0
        trace_count = 0
        while area_diff <= self.area_threshold:
            trace_count += 1
            fdir = self.flowdir_array[tuple(current_pixel)]
            if fdir == 0 or fdir == self.flowdir_ndv:
                break
            prev_pixel = current_pixel
            current_pixel = current_pixel + window[fdir]
            current_cell = self.identify_cell_by_pixel(current_pixel)
            if self.check_pixel_out_of_bounds(current_pixel) or self.check_cell_out_of_bounds(cell, current_cell):
                latest_outlet = prev_pixel
                break
            if self.cell_outlets[tuple(current_pixel)] == 1:
                latest_outlet = current_pixel
                current_area = self.flowacc_array[tuple(current_pixel)]
                area_diff = current_area - init_area
            elif self.check_cell_change(prev_pixel, current_pixel):
                latest_outlet = prev_pixel

        discharge_cell = self.identify_cell_by_pixel(latest_outlet)
        diff = discharge_cell - cell
        cell_dir = self.define_direction(diff)
        self.cells[tuple(cell)] = cell_dir

    def get_cell_outlet_flowacc(self, cell):
        cell_outlet = self.find_cell_outlet_by_mask(cell)
        outlet_flowacc = self.flowacc_array[tuple(cell_outlet)]
        return outlet_flowacc

    def fix_intersections(self):
        for i in range(self.mrows - 1):
            for j in range(self.mcols - 1):
                if self.cells[i + 1, j] == 128 and self.cells[i + 1, j + 1] == 32:
                    if self.get_cell_outlet_flowacc(np.array((i, j))) > self.get_cell_outlet_flowacc(
                            np.array((i, j + 1))):
                        self.cells[i + 1, j] = 64
                    else:
                        self.cells[i + 1, j + 1] = 64

                if self.cells[i, j] == 2 and self.cells[i, j + 1] == 8:
                    if self.get_cell_outlet_flowacc(np.array((i + 1, j + 1))) > self.get_cell_outlet_flowacc(
                            np.array((i + 1, j))):
                        self.cells[i, j + 1] = 4
                    else:
                        self.cells[i, j] = 4

                if self.cells[i, j + 1] == 8 and self.cells[i + 1, j + 1] == 32:
                    if self.get_cell_outlet_flowacc(np.array((i, j))) > self.get_cell_outlet_flowacc(
                            np.array((i + 1, j))):
                        self.cells[i, j + 1] = 16
                    else:
                        self.cells[i + 1, j + 1] = 16

                if self.cells[i, j] == 2 and self.cells[i + 1, j] == 128:
                    if self.get_cell_outlet_flowacc(np.array((i, j + 1))) > self.get_cell_outlet_flowacc(
                            np.array((i + 1, j + 1))):
                        self.cells[i, j] = 1
                    else:
                        self.cells[i + 1, j] = 1

    def execute(self):
        
        print('\tCreating null mask...')
        self.create_null_mask()
        
        print('\tAssigning outlet pixels...')
        self.assign_outlet_pixels()
        
        ncells = self.mrows * self.mcols
        cnt = 0
        interval = ncells // 100
        for i in range(self.mrows):
            for j in range(self.mcols):
                if cnt % interval == 0:
                    print(f'\t\tProcessing cell {cnt + 1} out of {ncells}')
                cnt += 1
                cell = np.array((i, j))
                if not self.check_null_cell(cell):
                    self.assign_cell_direction(cell)
        self.fix_intersections()