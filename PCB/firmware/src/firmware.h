#include <Arduino.h>
#include <Wire.h>
#include <DHT.h>
#include <RTClib.h>
#include <SPI.h>
#include <SD.h>
#include <avr/wdt.h>
#include <LiquidCrystal.h>

#include "sfm3019.h"

int sfmError = 0;
SfmConfig sfm4300;

#define TCAADDR 0x70

int sensorLoc[8] = {0, 0, 0, 0, 0, 0, 0, 0};

void tcaselect(uint8_t i) {
    if (i > 7) return;
    
    Wire.beginTransmission(TCAADDR);
    Wire.write(1 << i);
    Wire.endTransmission();  
}

void tcaUnSelect() {
    Wire.beginTransmission(TCAADDR);
    Wire.write(0);  // no channel selected
    Wire.endTransmission();
}

void scanI2C() {
        Serial.println("\nI2C Scanner");

    byte error, address;
    int nDevices;

    Serial.println("Scanning...");

    nDevices = 0;
    for(address = 1; address < 127; address++ ) {
        Wire.beginTransmission(address);
        error = Wire.endTransmission();

        if (error == 0) {
            Serial.print("I2C device found at address 0x");
            if (address<16) 
                Serial.print("0");
            Serial.print(address,HEX);
            Serial.println("  !");

            nDevices++;
        }
        else if (error==4) {
            Serial.print("Unknown error at address 0x");
            if (address<16) 
                Serial.print("0");
            Serial.println(address,HEX);
        }    
    }
    if (nDevices == 0)
        Serial.println("No I2C devices found\n");
    else
        Serial.println("done\n");

    // -------------

    Serial.println("\nTCA I2C Scanner");

for (uint8_t t=0; t<8; t++) {
    tcaselect(t);
    Serial.print("TCA Port #"); Serial.println(t);

    for (uint8_t addr = 0; addr<=127; addr++) {
    if (addr == TCAADDR) continue;

    Wire.beginTransmission(addr);
    if (!Wire.endTransmission()) {
        Serial.print("Found I2C 0x");  Serial.println(addr,HEX);
        if(addr == 0x31) {
            Serial.println("ELT Found!");
        }
        sensorLoc[t] = 1;
    }
    }
}
Serial.println("\ndone");

for(int i=0;i<8;i++) {
    Serial.println(sensorLoc[i]);
}


}



//// ELT Function (Y-Pod Method) ////

void wire_setup(int address, byte cmd, int from) 
{
	Wire.beginTransmission(address);
	Wire.write(cmd);
	Wire.endTransmission();
	Wire.requestFrom(address, from);
}

float getS300CO2()  
{
	int i = 1;
	long reading;
	//float CO2val;
    Wire.beginTransmission(0x31);
    Wire.write('E');
    Wire.endTransmission();
	wire_setup(0x31, 0x52, 7);

	while (Wire.available())
	{
	byte val = Wire.read();
	if (i == 2)  
	{
		reading = val;
		reading = reading << 8;
	}
	if (i == 3)  
	{
		reading = reading | val;
	}
	i = i + 1;
	}

	//Shift Calculation to Atheros
	//    CO2val = reading / 4095.0 * 5000.0;
	//    CO2val = reading;
	return reading;
}

/// From Sensirion Example Code //

SfmConfig initFlowSensor() {
  const char* driver_version = sfm_common_get_driver_version();
    if (driver_version) {
        // Serial.print("\nSFM driver version: ");
        // Serial.println(driver_version);
    } else {
        // Serial.println("fatal: Getting driver version failed");
    }

    /* Reset all I2C devices */
    sfmError = sensirion_i2c_general_call_reset();
    if (sfmError) {
        // Serial.println("General call reset failed");
    }

    /* Wait for the SFM3019 to initialize */
    sensirion_sleep_usec(SFM3019_SOFT_RESET_TIME_US);

    while (sfm3019_probe()) {
        Serial.println("SFM sensor probing failed");
        sensirion_sleep_usec(100000);
    }

    uint32_t product_number = 0;
    uint8_t serial_number[8] = {};
    sfmError = sfm_common_read_product_identifier(SFM3019_I2C_ADDRESS,
                                               &product_number, &serial_number);
    if (sfmError) {
        // Serial.println("Failed to read product identifier");
    } else {
        // Serial.print("product: 0x");
        // Serial.print(product_number, HEX);
        // Serial.print(" serial: 0x");
        for (size_t i = 0; i < 8; ++i) {
            // Serial.print(serial_number[i], HEX);
        }
        // Serial.println("");
    }
    
    sfm4300 = sfm3019_create();

    sfmError = sfm_common_start_continuous_measurement(
        &sfm4300, SFM3019_CMD_START_CONTINUOUS_MEASUREMENT_AIR);

    if (sfmError) {
        // Serial.println("Failed to start measurement");
        Serial.println("\t\toffline");
    } else {
        Serial.println("\t\tonline");
    }
 
    /* Wait for the first measurement to be available. Wait for
     * SFM3019_MEASUREMENT_WARM_UP_TIME_US instead for more reliable results */
    sensirion_sleep_usec(SFM3019_MEASUREMENT_INITIALIZATION_TIME_US);

    return sfm4300;
}

float readFlowSensor() {

  int16_t flow_raw;
  int16_t temperature_raw;
  uint16_t status;

  sfmError = sfm_common_read_measurement_raw(&sfm4300, &flow_raw,
                                          &temperature_raw, &status);
  if (sfmError) {
      Serial.println("Error while reading measurement.");
      return NULL;
  } else {
      float flow;
      float temperature;
      sfmError = sfm_common_convert_flow_float(&sfm4300, flow_raw, &flow);
      if (sfmError) {
          Serial.println("Error while converting flow");
          return NULL;
      }
      temperature = sfm_common_convert_temperature_float(temperature_raw);

      return flow;
  }
  
}
