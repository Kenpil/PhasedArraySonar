#include "mbed.h"
 
SPI spi(D11, D12, D13); // mosi, miso, sclk
DigitalOut cs(D10);//chip select?
Serial pc(USBTX, USBRX);


int main()
{
    
    spi.format(7, 0);
    spi.frequency(1000000);
 
    printf("Starting MCP3208\n");
    
    while(1)
    {
        cs = 0;
 
        spi.write(0x60);
        
        uint8_t d1 = spi.write(0x00);
        uint8_t d2 = spi.write(0x00);
        uint16_t d = (d1 << 8) | d2;
 
        float volt = 3.3 * (float)d / 4095;

        pc.printf("received\n");
        printf("d=%d temp=%5.2f\n", d, volt * 100);

        cs = 1;
        wait(1);
    }
}
