import matplotlib.pyplot as plt
import numpy as np

dados = np.load("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/mi_kd.npy")

complexos = list(dados[1:,0])
kd = list(dados[1:,2])
mi = list(dados[1:,1])

kd2 = []
for n in kd:
    n = float(n)
    n = round(n,3)
    kd2.append(n)
kd = kd2
mi2 = []
for n in mi:
    n = float(n)
    n = round(n,3)
    mi2.append(n)
mi = mi2

plt.scatter(mi, kd)
for pos, dt in enumerate(complexos):
    plt.annotate(dt, (mi[pos], kd[pos]))
plt.xlabel("Average MI")
plt.ylabel("Kd (uM)")

z = np.polyfit(mi, kd, 1)
p = np.poly1d(z)

plt.plot(mi, p(mi), "r-")

r = np.polyfit(mi, kd, 2)
rdois = np.poly1d(r)
yhat = rdois(mi)
media = np.sum(kd)/len(kd)
sst = np.sum((kd - media)**2)
ssreg = np.sum((yhat - media)**2)

r2 = ssreg / sst
plt.legend(["R2={:.5f}".format(r2)])
plt.show()  