# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import pylab as pl

# 快速傅立叶变换
sampling_rate = 100
fft_size = 802
x_axis = range(1, 803, 1) # 时间点
means_fourier = np.fft.rfft(means)/fft_size # 对实数信号进行FFT计算, 为了正确显示波形能量，还需要将rfft函数的结果除以fft_size
freqs = np.linspace(0, sampling_rate/2, fft_size/2+1) # rfft函数的返回值是N/2+1个复数，分别表示从0(Hz)到sampling_rate/2(Hz)的N/2+1点频率的成分, 频域上的x轴
# xfp = 20*np.log10(np.clip(np.abs(xf), 1e-20, 1e100)) # 计算每个频率分量的幅值，并通过 20*np.log10() 将其转换为以db单位的值。为了防止0幅值的成分造成log10无法计算，我们调用np.clip对xf的幅值进行上下限处理
means_fp = np.abs(means_fourier)


pl.figure(figsize=(8,4))
pl.subplot(211)
pl.plot(x_axis, means)
pl.xlabel(u"时间(秒)")
pl.title(u"156.25Hz和234.375Hz的波形和频谱")

pl.subplot(212)
pl.plot(freqs, means_fp)
pl.xlabel(u"频率(Hz)")
pl.subplots_adjust(hspace=0.4)
pl.show()