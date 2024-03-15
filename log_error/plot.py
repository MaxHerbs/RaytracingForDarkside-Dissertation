import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from scipy.optimize import curve_fit

f = open("output_variance.txt","r")

time_list = []
error_list = []

for line in f.readlines():
    data = line.replace("\n","").split(",")
    time_list.append(datetime.strptime(data[0], '%Y-%m-%d %H:%M:%S.%f'))
    error_list.append(float(data[1]))

time_list = np.array(time_list)
error_list = np.array(error_list)
start_time = time_list[0]
time_since_start = [(t - start_time).total_seconds() for t in time_list]
time_since_start = np.array(time_since_start)


time_since_start = time_since_start/(60*60)

time_since_start = time_since_start[4:]
error_list = error_list[4:]

plt.figure(figsize=(10, 6))

def fit(x, a, tau):
    return a * np.exp(-x / tau) 

bad_coords = np.where(time_since_start < 2)
fit_x = np.delete(time_since_start, bad_coords)
fit_y = np.delete(error_list, bad_coords)

param, cov = curve_fit(fit, fit_x, fit_y)
x_high_res = np.linspace(np.min(time_since_start), np.max(time_since_start),100)


plt.plot(x_high_res, fit(x_high_res, *param))
plt.plot(time_since_start, error_list)
plt.xlabel('Time since start (hours)')
plt.ylabel('Value')
plt.title('Scatter Plot of Time vs Value')
plt.xticks(rotation=45)
plt.tight_layout()
plt.clf()


def moving_average(data, window_size):
    cumsum = np.cumsum(data, dtype=float)
    cumsum[window_size:] = cumsum[window_size:] - cumsum[:-window_size]
    return cumsum[window_size - 1:] / window_size

smoothed_y = moving_average(error_list, window_size=5)
time_since_start_truncated = time_since_start[:len(smoothed_y)]

param, cov = curve_fit(fit, time_since_start_truncated, smoothed_y)
x_high_res = np.linspace(np.min(time_since_start), np.max(time_since_start),100)



def lin_fit(x,m,c):
    return m*x +c

log_param,log_cov = curve_fit(lin_fit,np.log(time_since_start_truncated),np.log(smoothed_y))
from math import sqrt

print("Log log parameters: {}".format(log_param))
trend_power = log_param[0]
trend_power_error = sqrt(log_cov[0][0])
print("Gradient = {:0.4f} +/- {:0.4f}".format(trend_power, trend_power_error))

def real_fit(x,a,c):
    return a * ((x)**(trend_power)) +c

real_param, cov = curve_fit(real_fit, time_since_start_truncated, smoothed_y)
print("Non-linear fit: {}".format(real_param))
print("A = {:0.3f} +/- {:0.3f}\nc = {:0.3f} +/- {:0.3f}".format(real_param[0], sqrt(cov[0][0]), real_param[1], sqrt(cov[1][1])))


plt.rcParams.update({'font.size': 28})
plt.grid()

plt.plot(np.log(time_since_start_truncated), lin_fit(np.log(time_since_start_truncated), *log_param))
plt.plot(np.log(time_since_start_truncated), np.log(smoothed_y))
plt.xlabel('Log(Duration of Simulation / hours)')
plt.ylabel('Log(Eror in Light\nVariation / %)')
plt.tight_layout()
plt.savefig("logerrorvsruntime.png")
plt.clf()

plt.grid()
plt.plot(time_since_start_truncated, real_fit(time_since_start_truncated, *real_param), linewidth=4, label="Polynomial Fit of Data")
plt.plot(time_since_start_truncated, smoothed_y, linewidth=3, label="Data")
plt.legend()
plt.xlabel('Duration of Simulation / hours')
plt.ylabel('Eror In Maximum\nLight Variation / %')
plt.tight_layout()
plt.savefig("errorvsruntime.png")
plt.clf()

# x_entries = -250
# plt.plot(time_since_start_truncated[x_entries:], real_fit(time_since_start_truncated[x_entries:], *real_param))
# plt.plot(time_since_start_truncated[x_entries:], smoothed_y[x_entries:])
# plt.show()

def calc_runtime(error):
    a = real_param[0]
    n = trend_power
    e = error
    c = real_param[1]
    one_perc_runtime = 10**(1/(-n) * np.log10(a/(e+c)))

    a_err = a + sqrt(cov[0][0])
    n_err = n +trend_power_error
    c_err = c + sqrt(cov[1][1])
    one_perc_runtime_err = (10**(1/(-n_err) * np.log10(a/(e+c))))**2 + (10**(1/(-n) * np.log10(a_err/(e+c))))**2 + (10**(1/(-n) * np.log10(a/(e+c_err))))**2
    one_perc_runtime_err = sqrt(one_perc_runtime_err) - one_perc_runtime

    return("Error: {}% --- Runtime: {} hours".format(error, "{:0.2f} +/- {:0.2f}".format(one_perc_runtime, np.abs(one_perc_runtime_err))))




print(calc_runtime(10))
print(calc_runtime(5))
print(calc_runtime(3))
print(calc_runtime(1))
print(calc_runtime(0.1))