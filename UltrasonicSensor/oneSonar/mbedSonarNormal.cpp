#include "mbed.h"

DigitalOut trigger(D6);
DigitalOut myled(LED1); //monitor trigger
DigitalOut myled2(LED2); //monitor echo
DigitalIn  echo(D7);
int ussdistance = 0;
int correction = 0;
Timer sonar;

int main()
{
    sonar.reset();
// measure actual software polling timer delays
// delay used later in time correction
// start timer
    sonar.start();
// min software polling delay to read echo pin
    while (echo==2) {};
    myled2 = 0;
// stop timer
    sonar.stop();
// read timer
    correction = sonar.read_us();
    printf("Approximate software overhead timer delay is %d uS\n\r",correction);
    //printf("ab%d \n\r",correction);
//Loop to read Sonar distance values, scale, and print
    while(1) {
// trigger sonar to send a ping
        trigger = 1;
        myled = 1;
        myled2 = 0;
        sonar.reset();
        {
            trigger = 1;
            wait_us(10.0);//default
            trigger = 0;
            wait_us(10.0);
        }
        myled = 0;
//wait for echo high
        while (echo==0) {};
        myled2=echo;
//echo high, so start timer
        sonar.start();
//wait for echo low
        while (echo==1) {};
//stop timer and read value
        sonar.stop();
//subtract software overhead timer delay and scale to cm
//distance = 10;
        ussdistance = (sonar.read_us()-correction)/58.0;
        myled2 = 0;
        printf(" %d cm \n\r",ussdistance);
//wait so that any echo(s) return before sending another ping
        wait(0.2);
    }
}
