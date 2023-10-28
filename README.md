# ECC - BCH code (15,7,2) BER simulation
以 C++ 撰寫的 BCH code 模擬程式
## simulation environment
1. SNR：0 : 0.5 : 12
2. 模擬次數：最少100,000次；根據100倍法則，總錯誤少於100個bits時不中斷。但節省模擬時間及BER=0的可能性，仍設最大次數設為3,000,000。
3. 接收訊號 $r[n]$ 以及雜訊 $w[n]$ 設定:
```math
r[n] = y[n] + w[n]
```
```math
w[n] ~ N(0,1) = \sqrt{\frac{N_0}{2 \times coderate}}
```
```math
N_0 = 10^\frac{SNR(dB)}{10}
```
## Explanation
紀錄於 Explanation.pdf

## Result
![This is an alt text.](/result.png "BER of BCH(15,7,2).")
