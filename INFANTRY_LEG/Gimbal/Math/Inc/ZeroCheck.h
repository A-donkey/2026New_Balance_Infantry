#ifndef __ZEROCHECK_H
#define __ZEROCHECK_H

typedef struct
{
	float Circle;	   //转过圈数
	float CountCycle;  //转过一圈的总计数周期
	float LastValue;   //检测过零量上一次的值
	float ActualValue; //检测过零量当前值
	float PreError;	   //检测量判断差值
} ZeroCheck_Typedef;

// 定义过零检测的方向模式
typedef enum {
    DIR_BOTH = 0,    // 双向检测（适用于云台、底盘等可随意正反转的机构）
    DIR_FORWARD,     // 仅正向检测（适用于只往一个方向打弹的拨盘、单向摩擦轮等）
    DIR_REVERSE      // 仅反向检测（根据你的电机安装方向决定）
} ZeroCheck_Mode_e;

float ZeroCheck(ZeroCheck_Typedef *Zero, float value, float CountCycle,ZeroCheck_Mode_e mode);

#endif
