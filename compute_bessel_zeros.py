from scipy.special import jn_zeros

with open("bessel_zeros.cfg", "w") as output:
  for i in range(11):
    zeros = jn_zeros(i, 50)
    for zero in zeros:
      output.write("{0}\t{1}\n".format(i, zero))
print("Finished!")
