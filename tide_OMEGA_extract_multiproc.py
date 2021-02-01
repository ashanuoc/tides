import glob
import numpy as np
from netCDF4 import Dataset
from scipy.io import savemat
from os import system, name
import os
import time as tic
import multiprocessing as mp


def wave_extract(year):

    files = glob.glob(filePath + 'f.e20.FXSD.f19_f19.001.cam.h1.' + str(year)+'[-0-9]*.nc')

    files.sort()

    # create nan arrays to allocate memory
    # geo-potential mean height

    gpz_mean = np.full([145, len(files)], np.nan)

    amp_OMEGA_DW1 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DW1 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DW2 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DW2 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DW3 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DW3 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DW4 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DW4 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DW5 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DW5 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DS0 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DS0 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DE1 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DE1 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DE2 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DE2 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DE3 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DE3 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DE4 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DE4 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_DE5 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_DE5 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_bg = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_bg = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SPW1 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SPW1 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SPW2 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SPW2 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SPW3 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SPW3 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SPW4 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SPW4 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SPW5 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SPW5 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SW1 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SW1 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SW2 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SW2 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SW3 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SW3 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SW4 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SW4 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SW5 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SW5 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SE1 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SE1 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SE2 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SE2 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SE3 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SE3 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SE4 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SE4 = np.full([nlat, nlevs, len(files)], np.nan)

    amp_OMEGA_SE5 = np.full([nlat, nlevs, len(files)], np.nan)
    phs_OMEGA_SE5 = np.full([nlat, nlevs, len(files)], np.nan)

    for days in range(len(files)):
        print('Processing {} day = {}\n'.format(year,days+1))

        nc_fid = Dataset(files[days], 'r', format="NETCDF4")
        time = nc_fid['time'][:].data
        lev = nc_fid['lev'][:].data
        lat = nc_fid['lat'][:].data

        # checking length of the time is 8. if not write to a text file
        if len(time) != 8:

            gpz_mean[:, days] = np.mean(nc_fid['Z3'][:].data, axis=(0,2,3))/1000

            print('skipped {} day = {}\n'.format(year,days+1))

            with open(dest + 'missing_time.txt', 'a') as f:
                w = files[days]+ '\t'+ 'time length: '+ str(len(time))+ '\n'
                f.write(w)
                continue

        OMEGA = nc_fid['OMEGA'][:].data

        gpz_mean[:,days] = np.mean(nc_fid['Z3'][:].data, axis=(0,2,3))/1000

        y = np.full([Nt, nlevs, nlat, Nz], np.nan, np.complex)
        # Pyy = np.full([Nt, nlevs, nlat, Nz], np.nan)
        ampl_s = np.full([Nt, nlevs, nlat, Nz], np.nan)
        ang_s = np.full([Nt, nlevs, nlat, Nz], np.nan)

        for i in range(nlat):
            for j in range(nlevs):

                y[:, j, i, :] = np.fft.fftshift(np.fft.fft2(OMEGA[:, j, i, :]))  # 2D fft
                ampl_s[:, j, i, :] = np.absolute(y[:, j, i, :])/(Nz*Nt)  # amplitude
                ang_s[:, j, i, :] = np.angle(y[:, j, i, :])  # phase

                # the amplitude of waves are doubled due to symmetry and only positive frequencies are considered

                amp_OMEGA_DW1[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == 1]
                phs_OMEGA_DW1[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == 1]

                amp_OMEGA_DW2[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == 2]
                phs_OMEGA_DW2[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == 2]

                amp_OMEGA_DW3[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == 3]
                phs_OMEGA_DW3[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == 3]

                amp_OMEGA_DW4[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == 4]
                phs_OMEGA_DW4[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == 4]

                amp_OMEGA_DW5[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == 5]
                phs_OMEGA_DW5[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == 5]

                amp_OMEGA_DS0[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == 0]
                phs_OMEGA_DS0[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == 0]

                amp_OMEGA_DE1[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == -1]
                phs_OMEGA_DE1[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == -1]

                amp_OMEGA_DE2[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == -2]
                phs_OMEGA_DE2[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == -2]

                amp_OMEGA_DE3[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == -3]
                phs_OMEGA_DE3[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == -3]

                amp_OMEGA_DE4[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == -4]
                phs_OMEGA_DE4[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == -4]

                amp_OMEGA_DE5[i, j, days] = 2 * ampl_s[freq == 1, j, i, wavenumber == -5]
                phs_OMEGA_DE5[i, j, days] = -ang_s[freq == 1, j, i, wavenumber == -5]

                amp_OMEGA_bg[i, j, days] = ampl_s[freq == 0, j, i, wavenumber == 0]
                phs_OMEGA_bg[i, j, days] = -ang_s[freq == 0, j, i, wavenumber == 0]

                amp_OMEGA_SPW1[i, j, days] = 2 * ampl_s[freq == 0, j, i, wavenumber == 1]
                phs_OMEGA_SPW1[i, j, days] = -ang_s[freq == 0, j, i, wavenumber == 1]

                amp_OMEGA_SPW2[i, j, days] = 2 * ampl_s[freq == 0, j, i, wavenumber == 2]
                phs_OMEGA_SPW2[i, j, days] = -ang_s[freq == 0, j, i, wavenumber == 2]

                amp_OMEGA_SPW3[i, j, days] = 2 * ampl_s[freq == 0, j, i, wavenumber == 3]
                phs_OMEGA_SPW3[i, j, days] = -ang_s[freq == 0, j, i, wavenumber == 3]

                amp_OMEGA_SPW4[i, j, days] = 2 * ampl_s[freq == 0, j, i, wavenumber == 4]
                phs_OMEGA_SPW4[i, j, days] = -ang_s[freq == 0, j, i, wavenumber == 4]

                amp_OMEGA_SPW5[i, j, days] = 2 * ampl_s[freq == 0, j, i, wavenumber == 5]
                phs_OMEGA_SPW5[i, j, days] = -ang_s[freq == 0, j, i, wavenumber == 5]

                amp_OMEGA_SW1[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == 1]
                phs_OMEGA_SW1[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == 1]

                amp_OMEGA_SW2[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == 2]
                phs_OMEGA_SW2[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == 2]

                amp_OMEGA_SW3[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == 3]
                phs_OMEGA_SW3[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == 3]

                amp_OMEGA_SW4[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == 4]
                phs_OMEGA_SW4[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == 4]

                amp_OMEGA_SW5[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == 5]
                phs_OMEGA_SW5[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == 5]

                amp_OMEGA_SE1[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == -1]
                phs_OMEGA_SE1[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == -1]

                amp_OMEGA_SE2[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == -2]
                phs_OMEGA_SE2[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == -2]

                amp_OMEGA_SE3[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == -3]
                phs_OMEGA_SE3[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == -3]

                amp_OMEGA_SE4[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == -4]
                phs_OMEGA_SE4[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == -4]

                amp_OMEGA_SE5[i, j, days] = 2 * ampl_s[freq == 2, j, i, wavenumber == -5]
                phs_OMEGA_SE5[i, j, days] = -ang_s[freq == 2, j, i, wavenumber == -5]
        
        

    DW1 = {"amp_OMEGA_DW1": amp_OMEGA_DW1, "phs_OMEGA_DW1": phs_OMEGA_DW1, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DW1/' + 'WACCMX_OMEGA_DW1_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DW1)

    DW2 = {"amp_OMEGA_DW2": amp_OMEGA_DW2, "phs_OMEGA_DW2": phs_OMEGA_DW2, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DW2/' + 'WACCMX_OMEGA_DW2_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DW2)

    DW3 = {"amp_OMEGA_DW3": amp_OMEGA_DW3, "phs_OMEGA_DW3": phs_OMEGA_DW3, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DW3/' + 'WACCMX_OMEGA_DW3_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DW3)

    DW4 = {"amp_OMEGA_DW4": amp_OMEGA_DW4, "phs_OMEGA_DW4": phs_OMEGA_DW4, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DW4/' + 'WACCMX_OMEGA_DW4_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DW4)

    DW5 = {"amp_OMEGA_DW5": amp_OMEGA_DW5, "phs_OMEGA_DW5": phs_OMEGA_DW5, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DW5/' + 'WACCMX_OMEGA_DW5_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DW5)

    DE1 = {"amp_OMEGA_DE1": amp_OMEGA_DE1, "phs_OMEGA_DE1": phs_OMEGA_DE1, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DE1/' + 'WACCMX_OMEGA_DE1_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DE1)

    DE2 = {"amp_OMEGA_DE2": amp_OMEGA_DE2, "phs_OMEGA_DE2": phs_OMEGA_DE2, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DE2/' + 'WACCMX_OMEGA_DE2_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DE2)

    DE3 = {"amp_OMEGA_DE3": amp_OMEGA_DE3, "phs_OMEGA_DE3": phs_OMEGA_DE3, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DE3/' + 'WACCMX_OMEGA_DE3_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DE3)

    DE4 = {"amp_OMEGA_DE4": amp_OMEGA_DE4, "phs_OMEGA_DE4": phs_OMEGA_DE4, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DE4/' + 'WACCMX_OMEGA_DE4_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DE4)

    DE5 = {"amp_OMEGA_DE5": amp_OMEGA_DE5, "phs_OMEGA_DE5": phs_OMEGA_DE5, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DE5/' + 'WACCMX_OMEGA_DE5_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DE5)

    DS0 = {"amp_OMEGA_DS0": amp_OMEGA_DS0, "phs_OMEGA_DS0": phs_OMEGA_DS0, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'DS0/' + 'WACCMX_OMEGA_DS0_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, DS0)

    SW1 = {"amp_OMEGA_SW1": amp_OMEGA_SW1, "phs_OMEGA_SW1": phs_OMEGA_SW1, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SW1/' + 'WACCMX_OMEGA_SW1_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SW1)

    SW2 = {"amp_OMEGA_SW2": amp_OMEGA_SW2, "phs_OMEGA_SW2": phs_OMEGA_SW2, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SW2/' + 'WACCMX_OMEGA_SW2_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SW2)

    SW3 = {"amp_OMEGA_SW3": amp_OMEGA_SW3, "phs_OMEGA_SW3": phs_OMEGA_SW3, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SW3/' + 'WACCMX_OMEGA_SW3_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SW3)

    SW4 = {"amp_OMEGA_SW4": amp_OMEGA_SW4, "phs_OMEGA_SW4": phs_OMEGA_SW4, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SW4/' + 'WACCMX_OMEGA_SW4_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SW4)

    SW5 = {"amp_OMEGA_SW5": amp_OMEGA_SW5, "phs_OMEGA_SW5": phs_OMEGA_SW5, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SW5/' + 'WACCMX_OMEGA_SW5_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SW5)

    SE1 = {"amp_OMEGA_SE1": amp_OMEGA_SE1, "phs_OMEGA_SE1": phs_OMEGA_SE1, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SE1/' + 'WACCMX_OMEGA_SE1_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SE1)

    SE2 = {"amp_OMEGA_SE2": amp_OMEGA_SE2, "phs_OMEGA_SE2": phs_OMEGA_SE2, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SE2/' + 'WACCMX_OMEGA_SE2_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SE2)

    SE3 = {"amp_OMEGA_SE3": amp_OMEGA_SE3, "phs_OMEGA_SE3": phs_OMEGA_SE3, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SE3/' + 'WACCMX_OMEGA_SE3_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SE3)

    SE4 = {"amp_OMEGA_SE4": amp_OMEGA_SE4, "phs_OMEGA_SE4": phs_OMEGA_SE4, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SE4/' + 'WACCMX_OMEGA_SE4_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SE4)

    SE5 = {"amp_OMEGA_SE5": amp_OMEGA_SE5, "phs_OMEGA_SE5": phs_OMEGA_SE5, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SE5/' + 'WACCMX_OMEGA_SE5_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SE5)

    SPW1 = {"amp_OMEGA_SPW1": amp_OMEGA_SPW1, "phs_OMEGA_SPW1": phs_OMEGA_SPW1, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SPW1/' + 'WACCMX_OMEGA_SPW1_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SPW1)

    SPW2 = {"amp_OMEGA_SPW2": amp_OMEGA_SPW2, "phs_OMEGA_SPW2": phs_OMEGA_SPW2, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SPW2/' + 'WACCMX_OMEGA_SPW2_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SPW2)

    SPW3 = {"amp_OMEGA_SPW3": amp_OMEGA_SPW3, "phs_OMEGA_SPW3": phs_OMEGA_SPW3, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SPW3/' + 'WACCMX_OMEGA_SPW3_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SPW3)

    SPW4 = {"amp_OMEGA_SPW4": amp_OMEGA_SPW4, "phs_OMEGA_SPW4": phs_OMEGA_SPW4, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SPW4/' + 'WACCMX_OMEGA_SPW4_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SPW4)

    SPW5 = {"amp_OMEGA_SPW5": amp_OMEGA_SPW5, "phs_OMEGA_SPW5": phs_OMEGA_SPW5, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'SPW5/' + 'WACCMX_OMEGA_SPW5_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, SPW5)

    OMEGA_bg = {"amp_OMEGA_bg": amp_OMEGA_bg, "phs_OMEGA_bg": phs_OMEGA_bg, "lev": lev, "lat": lat, "gpz_mean": gpz_mean}
    outfile_name = dest + 'OMEGA_bg/' + 'WACCMX_OMEGA_bg_short_term_one_day_' + str(year) + '.mat'
    savemat(outfile_name, OMEGA_bg)


if __name__ == "__main__":
    starttime = tic.time()

    filePath = '/data/avitharana/WACCMX/'  # origin
    dest = '/data/avitharana/WACCMX_Tides/Tides_OMEGA/'  # destination

    # Create folders if not exsist

    fld = ['DW1', 'DW2', 'DW3', 'DW4', 'DW5', 'DE1', 'DE2', 'DE3', 'DE4', 'DE5', 'DS0', 'SW1', 'SW2', 'SW3', 'SW4', 'SW5',
        'SE1', 'SE2', 'SE3', 'SE4', 'SE5', 'SPW1', 'SPW2', 'SPW3', 'SPW4', 'SPW5', 'OMEGA_bg']
    for item in fld:
        outpath = dest+item
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        else:
            print('already exist')

    # f = open(dest + 'missing_data.txt', "w")  # create file to save if a year of data not in folder
    # f.close()
    f = open(dest + 'missing_time.txt', "w")  # create file to save if time length is not 8
    f.close()

    Nt = 8  # Number of time steps in the 1-day window.

    nlat = 96
    nlevs = 145
    nsteps = 1  # number of days

    dz = 2.5  # Distance increment in degrees 360/64
    dt = 3  # time increment in Hours
    Nz = 144  # Number of samples available along longitude

    df = 1/(Nt*dt)  # temporal frequency
    dk = 1/(Nz*dz)  # spatial frequency

    wavenumber = (np.arange(1, Nz+1) - (Nz/2+1)) * dk * 360

    freq = (np.arange(1, Nt+1) - (Nt/2+1)) * df * 24  # we only use positive frequency here because of symmetry
    with np.errstate(divide='ignore'):
        freq_d = 1./freq

    num_workers = mp.cpu_count()
    pool = mp.Pool(num_workers)
    pool.map(wave_extract, range(1980, 2018))
    pool.close()
    endtime = tic.time()
    print(f"Time taken {endtime - starttime} seconds")