#include "LQR.h"

//定腿长下增益矩阵,仅用于判断极性
/*		
		{  -4.110908f,   -7.566761f,  -12.930175f,   -4.403529f,    8.670339f,    0.422445f,  -46.719406f,   -1.851758f,   82.142267f,    6.304709f},  // T_r_to_b
    {  -4.110908f,   -7.566761f,   12.930175f,    4.403529f,  -46.719406f,   -1.851758f,    8.670339f,    0.422445f,   82.142267f,    6.304709f},  // T_l_to_b
    {   1.819352f,    3.404594f,   -2.877683f,   -0.918506f,    2.887613f,    0.145660f,   19.865001f,    0.804548f,   17.550242f,    1.941559f},  // T_wr_to_r
    {   1.819352f,    3.404594f,    2.877683f,    0.918506f,   19.865001f,    0.804548f,    2.887613f,    0.145660f,   17.550242f,    1.941559f}   // T_wl_to_l
*/
//变腿长下多项式拟合的全部系数
float P[40][6] = 
{ 		
/*
lqr_Q = diag([200, 50, 1000, 50, 8000, 10, 8000, 10, 40000, 100]);
lqr_R = diag([2, 2, 20, 20]);
*/
		{    -8.98892f,      22.9346f,      4.13578f,     -14.2592f,     -38.0463f,       23.608f},  // K[0][0]
    {    -16.6244f,      45.5449f,      10.5408f,     -24.5723f,     -81.9709f,       41.412f},  // K[0][1]
    {    -8.87441f,     -23.0388f,     -29.9178f,       32.475f,     -14.3674f,      52.1503f},  // K[0][2]
    {    -1.96286f,     -5.52999f,     -13.4876f,      7.60431f,      -8.9573f,      21.4795f},  // K[0][3]
    {    -2.50401f,      39.2646f,      33.3319f,     -36.7014f,      19.4181f,     -45.2782f},  // K[0][4]
    {   -0.123389f,      1.68959f,      1.48838f,    -0.461997f,      2.64304f,      -2.6655f},  // K[0][5]
    {    -72.9834f,     -5.93523f,       157.92f,      18.1166f,     -40.5672f,     -148.242f},  // K[0][6]
    {    -2.24739f,    -0.349489f,      3.71373f,      1.07141f,      -2.6376f,     -4.59486f},  // K[0][7]
    {     48.4503f,     -103.988f,      350.395f,      135.953f,     -31.3588f,      -444.19f},  // K[0][8]
    {     2.95177f,     -7.27003f,      31.7343f,      11.1002f,     -10.7559f,     -34.3039f},  // K[0][9]
    {    -8.98892f,      4.13578f,      22.9346f,       23.608f,     -38.0463f,     -14.2592f},  // K[1][0]
    {    -16.6244f,      10.5408f,      45.5449f,       41.412f,     -81.9709f,     -24.5723f},  // K[1][1]
    {     8.87441f,      29.9178f,      23.0388f,     -52.1503f,      14.3674f,      -32.475f},  // K[1][2]
    {     1.96286f,      13.4876f,      5.52999f,     -21.4795f,       8.9573f,     -7.60431f},  // K[1][3]
    {    -72.9834f,       157.92f,     -5.93523f,     -148.242f,     -40.5672f,      18.1166f},  // K[1][4]
    {    -2.24739f,      3.71373f,    -0.349489f,     -4.59486f,      -2.6376f,      1.07141f},  // K[1][5]
    {    -2.50401f,      33.3319f,      39.2646f,     -45.2782f,      19.4181f,     -36.7014f},  // K[1][6]
    {   -0.123389f,      1.48838f,      1.68959f,      -2.6655f,      2.64304f,    -0.461997f},  // K[1][7]
    {     48.4503f,      350.395f,     -103.988f,      -444.19f,     -31.3588f,      135.953f},  // K[1][8]
    {     2.95177f,      31.7343f,     -7.27003f,     -34.3039f,     -10.7559f,      11.1002f},  // K[1][9]
    {     1.28076f,     -8.67226f,       16.456f,      5.04975f,      5.29176f,     -21.0383f},  // K[2][0]
    {      2.7783f,     -16.5136f,      26.2056f,      7.79677f,      15.9665f,     -37.7421f},  // K[2][1]
    {    -6.54469f,     -7.03601f,      25.7643f,      8.05418f,      9.55786f,     -29.5269f},  // K[2][2]
    {    -1.69147f,     -3.71824f,      7.82283f,      5.02094f,      1.55611f,     -7.04894f},  // K[2][3]
    {     2.25825f,      13.5467f,     -5.91166f,     -15.4044f,     -12.6186f,      5.74796f},  // K[2][4]
    {   0.0627841f,     0.836203f,    -0.292316f,    -0.442534f,      -1.4477f,     0.470086f},  // K[2][5]
    {     13.1939f,      -9.1897f,      79.4797f,      11.0404f,     0.628416f,     -90.4357f},  // K[2][6]
    {    0.395352f,    -0.407882f,      3.12074f,     0.596593f,    -0.164728f,     -1.83962f},  // K[2][7]
    {     33.3901f,     -34.6933f,     -66.3601f,      33.0457f,      32.9298f,      40.2323f},  // K[2][8]
    {     3.65134f,     -6.38399f,     -4.81626f,      5.75536f,      6.15308f,     0.483395f},  // K[2][9]
    {     1.28076f,       16.456f,     -8.67226f,     -21.0383f,      5.29176f,      5.04975f},  // K[3][0]
    {      2.7783f,      26.2056f,     -16.5136f,     -37.7421f,      15.9665f,      7.79677f},  // K[3][1]
    {     6.54469f,     -25.7643f,      7.03601f,      29.5269f,     -9.55786f,     -8.05418f},  // K[3][2]
    {     1.69147f,     -7.82283f,      3.71824f,      7.04894f,     -1.55611f,     -5.02094f},  // K[3][3]
    {     13.1939f,      79.4797f,      -9.1897f,     -90.4357f,     0.628416f,      11.0404f},  // K[3][4]
    {    0.395352f,      3.12074f,    -0.407882f,     -1.83962f,    -0.164728f,     0.596593f},  // K[3][5]
    {     2.25825f,     -5.91166f,      13.5467f,      5.74796f,     -12.6186f,     -15.4044f},  // K[3][6]
    {   0.0627841f,    -0.292316f,     0.836203f,     0.470086f,      -1.4477f,    -0.442534f},  // K[3][7]
    {     33.3901f,     -66.3601f,     -34.6933f,      40.2323f,      32.9298f,      33.0457f},  // K[3][8]
    {     3.65134f,     -4.81626f,     -6.38399f,     0.483395f,      6.15308f,      5.75536f}   // K[3][9]
};
float Offset_Fit_Coefficients[4] = {
					-17.178f,				11.939f,			-2.9224f,				0.3429f// theta_eq
};

//定义控制矩阵
float u[10];

//定义拟合增益矩阵
float Fitting_K[4][10];

void Offset_Calc(Flag_Bit_t *flag,Compensation_Amount_t *comp,float (*Fit_Coefficients),float TL){
 float theta_offset = Fit_Coefficients[0]*TL*TL*TL + Fit_Coefficients[1]*TL*TL +	 Fit_Coefficients[2]*TL	+	Fit_Coefficients[3];
 comp->Gravity_Comp_Theta_l = 	theta_offset;
 comp->Gravity_Comp_Theta_r = 	theta_offset;
 //theta_b定义正方向与开源不同
 if(flag->spinning_flag){
	comp->Spin_Comp_Theta_b = 0;
 }
 else{
	comp->Spin_Comp_Theta_b = 0;
 }
}

void Fitting_K_Calc(float (*fitting_k)[10],float (*p)[6],float L_l,float L_r)
{
	static unsigned short int i = 0;
	static unsigned short int j = 0;
	
	for(i=0;i<=3;i++)
	{
		for(j=0;j<=9;j++)
		{
			fitting_k[i][j] = p[i*10+j][0] + p[i*10+j][1]*L_l + p[i*10+j][2]*L_r + p[i*10+j][3]*(L_l*L_l) + p[i*10+j][4]*L_l*L_r + p[i*10+j][5]*(L_r*L_r);
		}
	}
}

float test_theta = 0.09f;
float test_yaw = 0.0f;
void LQR_Calc(Flag_Bit_t *flag,
							Goal_Setting_t *goal,
							Compensation_Amount_t *comp,
							Body_Current_Situation_t *body,
							Leg_Current_Situation_t *leg[2],
							Joint_Motor_Status_t *joint_m[2])
{
	//bz:|0,T_br|1,T_bl|2,T_wr|3,T_wl|
	static float T[4];
	u[0] = Max_Output(										   	 0  -        body->x,      X_MAX); u[1] = 	goal->d_x_t -       body->d_x;
	u[2] = Max_Output(	Find_Min_RADIAN(body->abs_yaw,goal->yaw_t),    Yaw_MAX); u[3] = goal->d_yaw_t -     body->d_yaw;
	u[4] = Max_Output(comp->Gravity_Comp_Theta_l  -  leg[0]->theta,Theta_L_MAX); u[5] =             0 - leg[0]->d_theta;
	u[6] = Max_Output(comp->Gravity_Comp_Theta_r  -  leg[1]->theta,Theta_R_MAX); u[7] =             0 - leg[1]->d_theta;
	u[8] = Max_Output(	 comp->Spin_Comp_Theta_b  -  	 body->theta,Theta_B_MAX); u[9] =             0 -   body->d_theta;
	
//	u[0] = Max_Output(										   	 0  -        body->x,      X_MAX); u[1] = 	goal->d_x_t -       body->d_x;
//	u[2] = Max_Output(	Find_Min_RADIAN(body->abs_yaw,goal->yaw_t-test_yaw * Ang_PI),    Yaw_MAX); u[3] = goal->d_yaw_t -     body->d_yaw;
//	u[4] = Max_Output(								test_theta  -  leg[0]->theta,Theta_L_MAX); u[5] =             0 - leg[0]->d_theta;
//	u[6] = Max_Output(								test_theta  -  leg[1]->theta,Theta_R_MAX); u[7] =             0 - leg[1]->d_theta;
//	u[8] = Max_Output(	 comp->Spin_Comp_Theta_b  -  	 body->theta,Theta_B_MAX); u[9] =             0 -   body->d_theta;
	
	//腿杆异常取消YAW增益提高腿杆THETA增益
	if(flag->theta_flag[0]||flag->theta_flag[1])
	{
		if(flag->theta_flag[0])
		{
			Fitting_K[0][2] = 0; Fitting_K[0][3] = 0; 					
			Fitting_K[2][2] = 0; Fitting_K[2][3] = 0; 	
			
			Fitting_K[0][4] = Fitting_K[0][4]*2.0f; Fitting_K[0][6] = Fitting_K[0][6]*2.0f; 
			Fitting_K[2][4] = Fitting_K[2][4]*2.0f; Fitting_K[2][6] = Fitting_K[2][6]*2.0f;
		}
		if(flag->theta_flag[1])
		{
			Fitting_K[1][2] = 0; Fitting_K[1][3] = 0; 					
			Fitting_K[3][2] = 0; Fitting_K[3][3] = 0; 	
			
			Fitting_K[1][4] = Fitting_K[1][4]*2.0f; Fitting_K[1][6] = Fitting_K[1][6]*2.0f; 
			Fitting_K[3][4] = Fitting_K[3][4]*2.0f; Fitting_K[3][6] = Fitting_K[3][6]*2.0f;
		}
	}

	//先调腿部摆角，再调机体PITCH，再调位移速度，最后调YAW
		if(!flag->off_flag) 
		{
			if(!flag->spinning_flag) 
			{
					//关节电机极性
				T[0] = 
				-u[0]*Fitting_K[0][0]-u[1]*Fitting_K[0][1]
				+u[2]*Fitting_K[0][2]-u[3]*Fitting_K[0][3]
				-u[4]*Fitting_K[0][4]-u[5]*Fitting_K[0][5]
				-u[6]*Fitting_K[0][6]-u[7]*Fitting_K[0][7]
				-u[8]*Fitting_K[0][8]-u[9]*Fitting_K[0][9];

				T[1] = 
				-u[0]*Fitting_K[1][0]-u[1]*Fitting_K[1][1]
				+u[2]*Fitting_K[1][2]-u[3]*Fitting_K[1][3]
				-u[4]*Fitting_K[1][4]-u[5]*Fitting_K[1][5]
				-u[6]*Fitting_K[1][6]-u[7]*Fitting_K[1][7]
				-u[8]*Fitting_K[1][8]-u[9]*Fitting_K[1][9];
				
					//轮电机极性
				T[2] = 
				+u[0]*Fitting_K[2][0]+u[1]*Fitting_K[2][1]
				-u[2]*Fitting_K[2][2]+u[3]*Fitting_K[2][3]
				+u[4]*Fitting_K[2][4]+u[5]*Fitting_K[2][5]
				+u[6]*Fitting_K[2][6]+u[7]*Fitting_K[2][7]
				+u[8]*Fitting_K[2][8]+u[9]*Fitting_K[2][9];
					
				T[3] = 
				+u[0]*Fitting_K[3][0]+u[1]*Fitting_K[3][1]
				-u[2]*Fitting_K[3][2]+u[3]*Fitting_K[3][3]
				+u[4]*Fitting_K[3][4]+u[5]*Fitting_K[3][5]
				+u[6]*Fitting_K[3][6]+u[7]*Fitting_K[3][7]
				+u[8]*Fitting_K[3][8]+u[9]*Fitting_K[3][9];
				
			}
			else //小陀螺取消X位移,YAW位置增益
			{
				
					//关节电机极性
				T[0] = 
				- 0												 - u[1]*Fitting_K[0][1]
				- 0												 - u[3]*Fitting_K[0][3]
				- u[4]*Fitting_K[0][4]*1.5f- u[5]*Fitting_K[0][5]*1.5f
				-	u[6]*Fitting_K[0][6]*1.5f- u[7]*Fitting_K[0][7]*1.5f
				- u[8]*Fitting_K[0][8]*2.0f- u[9]*Fitting_K[0][9]*2.0f;

				T[1] = 
				- 0												 - u[1]*Fitting_K[1][1]
				- 0												 - u[3]*Fitting_K[1][3]
				- u[4]*Fitting_K[1][4]*1.5f- u[5]*Fitting_K[1][5]*1.5f
				-	u[6]*Fitting_K[1][6]*1.5f- u[7]*Fitting_K[1][7]*1.5f
				- u[8]*Fitting_K[1][8]*2.0f- u[9]*Fitting_K[1][9]*2.0f;
				
					//轮电机极性
				T[2] = 
				+	0												 + u[1]*Fitting_K[2][1]
				+ 0												 + u[3]*Fitting_K[2][3]
				+ u[4]*Fitting_K[2][4]*1.5f+ u[5]*Fitting_K[2][5]*1.5f
				+	u[6]*Fitting_K[2][6]*1.5f+ u[7]*Fitting_K[2][7]*1.5f
				+ u[8]*Fitting_K[2][8]*2.0f+ u[9]*Fitting_K[2][9]*2.0f;

				T[3] = 
				+	0												 + u[1]*Fitting_K[3][1]
				+ 0												 + u[3]*Fitting_K[3][3]
				+ u[4]*Fitting_K[3][4]*1.5f+ u[5]*Fitting_K[3][5]*1.5f
				+	u[6]*Fitting_K[3][6]*1.5f+ u[7]*Fitting_K[3][7]*1.5f
				+ u[8]*Fitting_K[3][8]*2.0f+ u[9]*Fitting_K[3][9]*2.0f;
				
			}
		}
		else //离地时只控制腿杆垂直
		{
			
					//关节电机极性
				T[0] = 
				- u[4]*Fitting_K[0][4]- u[5]*Fitting_K[0][5]
				-	u[6]*Fitting_K[0][6]- u[7]*Fitting_K[0][7];

				T[1] = 
				- u[4]*Fitting_K[1][4]- u[5]*Fitting_K[1][5]
				-	u[6]*Fitting_K[1][6]- u[7]*Fitting_K[1][7];

					//轮电机极性
				T[2] = 
				+ u[4]*Fitting_K[2][4]+ u[5]*Fitting_K[2][5]
				+	u[6]*Fitting_K[2][6]+ u[7]*Fitting_K[2][7];

				T[3] = 
				+ u[4]*Fitting_K[3][4]+ u[5]*Fitting_K[3][5]
				+	u[6]*Fitting_K[3][6]+ u[7]*Fitting_K[3][7];
			
		}

	//输出赋值
	joint_m[0]->Tp 							= T[1];joint_m[1]->Tp 							= T[0]; 
	joint_m[0]->T_Wheel     		= T[3];joint_m[1]->T_Wheel      		= T[2];
}








