# ECC - BCH code (15,7,2) BER simulation
以 C++ 撰寫的 BCH code 模擬程式
## simulation environment
1. SNR：0 : 0.5 : 12
2. 模擬次數：最少100,000次；根據100倍法則，總錯誤少於100個bits時不中斷。但節省模擬時間及BER=0的可能性，仍設最大次數設為3,000,000。
3. 接收訊號$ 行內公式 $ r[n] = y[n] + w[n],  w[n] is noise $ ，0 ~N(0,1) = 2*Nnoisecoderate，__10010SNRindBN
## Function
1. RIS_controller.v : main file
2. _7Seg.v          : 七段顯示器的顯示電路 (7-segment)
3. fDIV_115200HZ.v  : 除頻器 (更改 clk 頻率至 115200Hz)
4. one_shot.v       : 防止電路訊號跳動，造成錯誤的正負緣觸發
5. ROM_GPIO.v       : RIS 控制狀態編/解碼
6. SNRcollector.v   : SNR / 控制狀態 (serial-in data) 儲存、parallel-out
7. TxEncoder.v      : 將欲回傳資訊編碼
8. Tx.v             : RS232 Transmitter (115200Hz)
9. Rx.v             : RS232 Receiver (115200Hz)
## System Block Diagram
**Whole system block generate by RTL viewer**
![This is an alt text.](/Image/Whole_System.png)
**Rx system block generate by RTL viewer**
![This is an alt text.](/Image/Rx.png)
**其餘自動產生的細節可參考/Image/Whole_System(detail).pdf 圖可放大**
