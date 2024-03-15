import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


pos, var, err = [], [], []

f = open("output_variance.txt", "r")
for line in f.readlines():
    data = line.replace("\n","").split(",")
    pos.append(float(data[0]))
    var.append(float(data[1]))
    err.append(float(data[2]))


pos = np.array(pos)
var = np.array(var)
err = np.array(err)
pos = pos - 0.22
pos = 100*pos



negative_displacement_pos = pos[np.where(pos <= 0)]
negative_displacement_var = var[np.where(pos <= 0)]
negative_displacement_err = err[np.where(pos <= 0)]

positive_displacement_pos = pos[np.where(pos >= 0)]
positive_displacement_var = var[np.where(pos >= 0)]
positive_displacement_err = err[np.where(pos >= 0)]



def fit_line(x, m,c):
    return m*x+c


def cubic(x, a,b,c,d):
    return a*x**3 + b*x**2 +c*x + d

plt.rcParams.update({'font.size': 22})
plt.subplots(2, 1, gridspec_kw={'height_ratios': [5, 1]}, figsize=(8, 10))
plt.subplot(2,1,1)
negative_fit, negative_cov = curve_fit(fit_line, negative_displacement_pos, negative_displacement_var, sigma=negative_displacement_err)
plt.plot(np.sort(negative_displacement_pos), fit_line(np.sort(negative_displacement_pos), *negative_fit), label="Linear fit of negative positions")


positive_fit, positive_cov = curve_fit(fit_line, positive_displacement_pos, positive_displacement_var, sigma=positive_displacement_err)
plt.plot(np.sort(positive_displacement_pos), fit_line(np.sort(positive_displacement_pos), *positive_fit), label="Linear fit of positive positions")


cubic_param, cov = curve_fit(cubic, pos,var,sigma=err)
x = np.linspace(np.min(pos), np.max(pos), 100)
#plt.plot(x, cubic(x,*cubic_param))


def compount_fit(x,y,neg_param,pos_param):
    x_return = []
    y_return = []
    for i in range(len(x)):
        if x[i] < 0 :
            x_return.append(x[i])
            y_return.append(neg_param[0] * x[i] + neg_param[1])
        if x[i] >=0:
            x_return.append(x[i])
            y_return.append(pos_param[0] * x[i] + pos_param[1])
    return np.array(x_return), np.array(y_return)



comp_fit = compount_fit(pos, var, negative_fit, positive_fit)
print(pos)
print(var)
print(err)
print(comp_fit)

resid = (var-comp_fit[1])/(err)





plt.errorbar(pos,var, err, fmt="o", label="Data", capsize=3)
plt.legend(fontsize=18)
plt.ylabel("Percentage variation from\nmaximum relative power")

plt.subplot(2, 1, 2)  
plt.scatter(pos, resid)
plt.ylim((-1,1))
plt.xlabel("Distance from top of PDU/cm")
plt.ylabel("Norm.\nResiduals")
plt.axhline(y=0, color='black', linestyle='--', linewidth=1)

plt.savefig("power_variation_with_position.png",  bbox_inches='tight')

from math import sqrt
output_txt = "### FIT DATA ###\n\nNegative fit\nm = {:0.2f}±{:0.2f}\nc = {:0.2f}±{:0.2f}\n\nPositive fit\nm = {:0.2f}±{:0.2f}\nc = {:0.2f}±{:0.2f}".format(negative_fit[0], sqrt(negative_cov[0][0]),negative_fit[1], sqrt(negative_cov[1][1]),positive_fit[0], sqrt(positive_cov[0][0]),positive_fit[1], sqrt(positive_cov[1][1]))
print(output_txt)

m1 = negative_fit[0]
c1 = negative_fit[1]
m2 = positive_fit[0]
c2 = positive_fit[1]



intersect = ((c2-c1)/m2)/(m1/m2+1)

from math import sqrt
m1_error  = m1 + sqrt(negative_cov[0][0])
c1_error = c1 + sqrt(negative_cov[1][1])
m2_error = m2 + sqrt(positive_cov[0][0])
c2_error = c2 + sqrt(positive_cov[1][1])

intersect_err = (((c2_error-c1)/m2)/(m1/m2+1))**2 + (((c2-c1_error)/m2)/(m1/m2+1))**2 + (((c2-c1)/m2)/(m1_error/m2+1))**2 + (((c2-c1)/m2)/(m1/m2_error+1))**2
intersect_err = sqrt(intersect_err) - intersect
print("Intersect = {} +/- {}".format(intersect, intersect_err))

variance = m1 * intersect + c1
variance_err = (m1_error * intersect + c1)**2 + (m1 * intersect + c1_error)**2
variance_err = sqrt(variance_err) - variance

print(f"Variance: {variance} +/- {variance_err}")