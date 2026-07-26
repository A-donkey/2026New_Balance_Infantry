#include "Super_Cap.h"
#include "Power_Limit.h"

//全局变量定义部分
volatile pm_od_t pm_od;

/*******************************************************************************************************
设定功率值
********************************************************************************************************/
void Pm_Power_Set(CAN_HandleTypeDef *hcan,uint16_t power,uint16_t buffer_energy)
{
	int16_t power_predict = P_cmd * 100;
	uint8_t data[8];
	
	data[0] = (uint8_t)(buffer_energy &  0xFF);
	data[1] = (uint8_t)(buffer_energy >>    8);
	data[2] = (uint8_t)(power &  0xFF);
	data[3] = (uint8_t)(power >>    8);
	data[4] = (uint8_t)(power_predict &  0xFF);
	data[5] = (uint8_t)(power_predict >>    8);
	data[6] = 0;
	data[7] = 0;

	Can_TxMessage(hcan,0x50,0x08,data);
}
