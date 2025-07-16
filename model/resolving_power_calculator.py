import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from astropy.io import fits
from scipy import signal
import sys
sys.path.append("/home/sborthakur/Documents/PhD_Tartu1/pyreduce_learning/PyReduce")
from pyreduce.util import polyfit2d
from scipy.optimize import curve_fit


def gen_image(data):
    vmin_init = np.mean(data) - np.std(data)
    vmax_init = np.mean(data) + np.std(data)

    # Create the figure and axis
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)  # space for slider

    # Display the image
    im = ax.imshow(data, cmap='gray', vmin=vmin_init, vmax=vmax_init, aspect="auto")
    plt.colorbar(im, ax=ax)

    # Add sliders for vmin and vmax (contrast control)
    ax_vmin = plt.axes([0.2, 0.1, 0.65, 0.03])
    ax_vmax = plt.axes([0.2, 0.05, 0.65, 0.03])
    slider_vmin = Slider(ax_vmin, 'Min', np.mean(data) - 5 * np.std(data),
                         np.mean(data), valinit=vmin_init)
    slider_vmax = Slider(ax_vmax, 'Max', np.mean(data),
                         np.mean(data) + 5 * np.std(data), valinit=vmax_init)

    # Update function
    def update(val):
        im.set_clim(vmin=slider_vmin.val, vmax=slider_vmax.val)
        fig.canvas.draw_idle()

    # Connect sliders to update function
    slider_vmin.on_changed(update)
    slider_vmax.on_changed(update)

    plt.show()

def create_image_from_lines(lines,ncol):
    min_order = int(np.min(lines["order"]))
    max_order = int(np.max(lines["order"]))

    img = np.zeros((max_order + 1, ncol))
    for line in lines:
        if line["order"] < 0:
            continue
        if line["xlast"] < 0 or line["xfirst"] > ncol:
            continue
        first = int(max(line["xfirst"], 0))
        last = int(min(line["xlast"], ncol))
        img[int(line["order"]), first:last] = line[
            "height"
        ] * signal.windows.gaussian(last - first, line["width"])
    return img

def apply_alignment_offset(lines, offset, select=None):
    if select is None:
        select = slice(None)
    lines["xfirst"][select] += offset[1]
    lines["xlast"][select] += offset[1]
    lines["posm"][select] += offset[1]
    lines["order"][select] += offset[0]
    return lines

def correlate_linelist(data,lines):
    nord, ncol = data.shape
    img = create_image_from_lines(lines, ncol)
    data = (data - np.mean(data)) / np.std(data)
    img = (img - np.mean(img)) / np.std(img)
    correlation = signal.correlate2d(data, img, mode="same")
    offset_order, offset_x = np.unravel_index(np.argmax(correlation), correlation.shape)
    # offset_order,offset_x = 22,1064
    offset_order = offset_order - img.shape[0] // 2 + 1 + 1
    offset_x = offset_x - img.shape[1] // 2 + 1
    offset = [int(offset_order), int(offset_x)]
    # apply offset
    lines = apply_alignment_offset(lines, offset)
    return lines
##############################################################################
def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-0.5 * ((x - mean) / stddev) ** 2)


def gen_data_from_npz(params, shape):
    data = np.zeros((max(params['order']), shape))
    ords = params['order'] - 1
    xfirst = params['xfirst']
    xlast = params['xlast']
    for i in range(ords.shape[0]):
        data[ords[i], xfirst[i]:xlast[i] + 1] = 1
    return ords, xfirst, xlast, data

def gen_wave_2d_image(linelist, deg, nord, ncol):
    params = linelist['cs_lines']
    x = params['posm']
    y = params['order']
    z = params['wlc']

    coeff = polyfit2d(x, y, z, deg)

    y_grid, x_grid = np.indices((nord, ncol))
    wave_img = np.polynomial.polynomial.polyval2d(x_grid, y_grid, coeff)
    return wave_img


def sigma_clipping(x, y, n_iter=20, sigma_lower=1, sigma_upper=3, plot=True, **kwargs):
    x_orig, y_orig = x, y
    if 'fontsize' in kwargs.keys():
        fontsize = kwargs['fontsize']
    else:
        fontsize = 20
    if plot:
        fig = plt.figure(figsize=(10, 10))
    for i in range(n_iter):
        yfit, _ = gaussian_fit(x, y)
        ynorm = y / yfit
        ymean = np.median(ynorm)
        ystd = np.std(ynorm)
        ind = np.where((ynorm >= ymean - ystd * sigma_lower) & (ynorm <= ymean + ystd * sigma_upper))[0]
        x, y = x[ind], y[ind]

        if plot:
            plt.clf()
            plt.plot(x_orig, y_orig, 'k.', alpha=0.7, markersize=0.5)
            plt.plot(x, y, 'r.', alpha=0.7, markersize=0.5)
            plt.title(f"Iteration - {i}", fontsize=fontsize)
            plt.xlabel("Velocity(km/s)", fontsize=fontsize)
            plt.ylabel("Arbitrary Flux", fontsize=fontsize)
            plt.xticks(fontsize=fontsize)
            plt.yticks(fontsize=fontsize)
            plt.tight_layout()
            plt.draw()
            plt.pause(0.5)
    if plot:
        plt.xlabel("Velocity(km/s)", fontsize=fontsize)
        plt.ylabel("Arbitrary Flux", fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.tight_layout()
        plt.show()

    return x, y

def makecombline(spec2d, wave2d, xfirst, xlast, ords, **kwargs):
    min_ord, max_ord = 0, spec2d.shape[0]
    min_wave_ind, max_wave_ind = 0, wave2d.shape[1]
    arraylen = 200

    if 'ord_range' in kwargs.keys():
        min_ord, max_ord = kwargs['ord_range'][0], kwargs['ord_range'][1]

    if 'wave_range_ind' in kwargs.keys():
        min_wave_ind, max_wave_ind = kwargs['wave_range_ind'][0], kwargs['wave_range_ind'][1]

    if "plotobject" in kwargs.keys():
        axs = kwargs["plotobject"]

    if "arraylen" in kwargs.keys():
        arraylen = int(kwargs["arraylen"])
    spec2d = spec2d[min_ord:max_ord+1,min_wave_ind:max_wave_ind+1]
    wave2d = wave2d[min_ord:max_ord+1,min_wave_ind:max_wave_ind+1]

    ind = np.where((ords>=min_ord)&(ords<max_ord+1)&(xfirst>=min_wave_ind)&(xfirst<max_wave_ind)&(xlast>=min_wave_ind)&(xlast<max_wave_ind))[0]
    ords,xfirst,xlast = ords[ind],xfirst[ind],xlast[ind]

    ords,xfirst,xlast = ords-min_ord, xfirst-min_wave_ind, xlast-min_wave_ind
    maxvel = 30
    velrange = np.linspace(-maxvel, maxvel, arraylen)
    combspec = np.zeros(len(velrange))
    sumweights = []

    num_lines = 0
    for ord,xi,xf in zip(ords,xfirst,xlast):
        speccrop,wavecrop = spec2d[ord,xi:xf+1],wave2d[ord,xi:xf+1]
        imax = np.argmax(speccrop)
        wmin = wavecrop[0]
        wmax = wavecrop[-1]
        initial_guess = [np.mean(speccrop), wavecrop[imax],
                         (wmax - wmin) / 2]

        params, _ = curve_fit(gaussian, wavecrop, speccrop, p0=initial_guess)
        amp, mean, stddev = params

        speccropfit = gaussian(wavecrop, amp, mean, stddev)

        c = 3e5

        vel, noise = c * (wavecrop - mean) / wavecrop, ((speccrop - speccropfit) / speccropfit) ** 2

        initial_guess = [np.mean(speccrop), 0, 2]
        params, _ = curve_fit(gaussian, vel, speccrop, p0=initial_guess)
        amp, mean, stddev = params
        resolving_power = int(0.5 * c / np.abs(stddev))
        weight = 1 / np.std(speccrop / speccropfit)
        if weight<0.03 and resolving_power>20_000:
            num_lines+=1
            combspec += weight * gaussian(velrange, amp, mean, stddev) / gaussian(mean, amp, mean, stddev)
            sumweights.append(weight)

    sumweights = np.array(sumweights)
    comblinesvel, combllinesnormflux = velrange, combspec / np.sum(sumweights)

    return comblinesvel, combllinesnormflux, num_lines

def gaussian_fit(x, y):
    xmax = x[-1]
    xmin = x[0]
    initial_guess = [np.mean(y), 0, (xmax - xmin) / 2]
    bounds = ([-np.inf, -np.inf, 0], [np.inf, np.inf, np.inf])
    params, _ = curve_fit(gaussian, x, y, p0=initial_guess, bounds=bounds)
    amp, mean, stddev = params

    yfit = gaussian(x, amp, mean, stddev)
    return yfit, params
###################################################################################

def resolving_power_calculator(thar_master_path,linelist_path,ndim,deg=[5,5],arraylen = 200):

    # ThAr lamp flux in pixels data
    hdu = fits.open(thar_master_path)
    thar_master_data = hdu[0].data
    nord,ncol = thar_master_data.shape

    # ThAr linelist data
    linelist = np.load(linelist_path,allow_pickle=True)
    linelist_params = linelist['cs_lines']
    linelist_params = correlate_linelist(thar_master_data,linelist_params)


    shape = ncol
    ords,xfirst,xlast,linelist_data = gen_data_from_npz(linelist_params,shape)
    wave2d = gen_wave_2d_image(linelist,deg,nord,ncol)

    nord_arr = np.arange(0, ndim + 1) * int((nord - 1) / ndim)
    ncol_arr = np.arange(0, ndim + 1) * int((ncol - 1) / ndim)

    nrows = len(nord_arr) - 1
    ncols = len(ncol_arr) - 1
    fig, axs = plt.subplots(nrows, ncols, figsize=(15, 12), sharex=True, sharey=True)
    axs = np.atleast_2d(axs)

    fullvel = np.zeros((nrows,ncols,arraylen))
    fullnormflux = np.zeros((nrows,ncols,arraylen))
    fullnormfluxfit = np.zeros((nrows,ncols,arraylen))
    fullresolving_power = np.zeros((nrows,ncols,1))

    for i in range(nrows):
        for j in range(ncols):
            comblinesvel, combllinesnormflux, num_lines = makecombline(
                thar_master_data,
                wave2d,
                xfirst,
                xlast,
                ords,
                ord_range=[nord_arr[i], nord_arr[i + 1]],
                wave_range_ind=[ncol_arr[j], ncol_arr[j + 1]],
                arraylen = arraylen
            )

            combllinesnormfluxfit, params = gaussian_fit(comblinesvel, combllinesnormflux)
            amp, mean, stddev = params

            c = 3e5
            resolving_power = int(0.5 * c / np.abs(stddev))
            fullvel[i][j],fullnormflux[i][j],fullnormfluxfit[i][j],fullresolving_power[i][j] = comblinesvel,combllinesnormflux,combllinesnormfluxfit,resolving_power

    return fullvel,fullnormflux,fullnormfluxfit,fullresolving_power