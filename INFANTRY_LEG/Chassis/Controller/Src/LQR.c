#include "LQR.h"

//定腿长下增益矩阵,仅用于判断极性
/*		
		{  -2.953531f,   -4.160691f,  -26.092797f,   -4.985259f,    8.084860f,    0.391328f,  -25.599030f,   -1.717345f,   57.058423f,    4.094747f},  // T_r_to_b
    {  -2.953531f,   -4.160691f,   26.092797f,    4.985259f,  -25.599030f,   -1.717345f,    8.084860f,    0.391328f,   57.058423f,    4.094747f},  // T_l_to_b
    {   1.275800f,    1.847248f,   -5.649478f,   -0.775354f,    2.646254f,    0.118954f,   10.072200f,    0.693241f,   12.740671f,    1.424436f},  // T_wr_to_r
    {   1.275800f,    1.847248f,    5.649478f,    0.775354f,   10.072200f,    0.693241f,    2.646254f,    0.118954f,   12.740671f,    1.424436f}   // T_wl_to_l
*/
//变腿长下多项式拟合的全部系数
float P[40][6] = 
{ 		
/*
lqr_Q = diag([100, 1, 4000, 1, 1000, 10, 1000, 10, 20000, 1]);
lqr_R = diag([2, 2, 20, 20]);
*/

		{    -5.27476f,      9.16828f,      5.98917f,     -5.80975f,     -13.5511f,      4.37555f},  // K[0][0]
    {     -7.6063f,      14.4019f,      9.45927f,     -8.79592f,     -21.5236f,       4.8074f},  // K[0][1]
    {    -13.8599f,     -38.7543f,     -53.3359f,      47.2595f,     -19.6461f,      95.5732f},  // K[0][2]
    {    -1.60501f,     -9.38614f,     -15.0846f,      10.7926f,     -5.51744f,      24.3972f},  // K[0][3]
    {    -5.18549f,      44.1868f,      44.1111f,       -33.79f,      26.7799f,     -62.7207f},  // K[0][4]
    {   -0.189302f,      1.59574f,      1.59505f,     0.507618f,      2.60913f,     -2.90095f},  // K[0][5]
    {    -25.1257f,      -12.854f,      2.19213f,       24.463f,     -48.3074f,      24.7253f},  // K[0][6]
    {    -2.08855f,    -0.484016f,      3.65408f,      1.25538f,     -3.03029f,     -4.80353f},  // K[0][7]
    {     31.8408f,     -64.5997f,      242.136f,      77.1085f,     -24.6292f,     -286.866f},  // K[0][8]
    {     1.96458f,      -7.0368f,       22.134f,      8.28557f,     -4.66542f,     -23.4424f},  // K[0][9]
    {    -5.27476f,      5.98917f,      9.16828f,      4.37555f,     -13.5511f,     -5.80975f},  // K[1][0]
    {     -7.6063f,      9.45927f,      14.4019f,       4.8074f,     -21.5236f,     -8.79592f},  // K[1][1]
    {     13.8599f,      53.3359f,      38.7543f,     -95.5732f,      19.6461f,     -47.2595f},  // K[1][2]
    {     1.60501f,      15.0846f,      9.38614f,     -24.3972f,      5.51744f,     -10.7926f},  // K[1][3]
    {    -25.1257f,      2.19213f,      -12.854f,      24.7253f,     -48.3074f,       24.463f},  // K[1][4]
    {    -2.08855f,      3.65408f,    -0.484016f,     -4.80353f,     -3.03029f,      1.25538f},  // K[1][5]
    {    -5.18549f,      44.1111f,      44.1868f,     -62.7207f,      26.7799f,       -33.79f},  // K[1][6]
    {   -0.189302f,      1.59505f,      1.59574f,     -2.90095f,      2.60913f,     0.507618f},  // K[1][7]
    {     31.8408f,      242.136f,     -64.5997f,     -286.866f,     -24.6292f,      77.1085f},  // K[1][8]
    {     1.96458f,       22.134f,      -7.0368f,     -23.4424f,     -4.66542f,      8.28557f},  // K[1][9]
    {    0.681807f,     -3.98993f,      8.17753f,      3.71218f,     0.947507f,     -10.1914f},  // K[2][0]
    {     1.13808f,     -5.94096f,      10.5199f,      5.39406f,      2.08977f,     -13.6433f},  // K[2][1]
    {    -10.5617f,     -10.3498f,      41.9223f,      10.2831f,      22.4015f,     -51.5766f},  // K[2][2]
    {    -1.74492f,     -2.42508f,      8.30697f,      2.25218f,      5.43811f,     -9.30878f},  // K[2][3]
    {     2.70925f,      12.1661f,     -7.42684f,     -9.89115f,     -26.4818f,      9.94957f},  // K[2][4]
    {   0.0435913f,     0.888189f,    -0.245292f,    -0.453017f,     -2.04902f,     0.625218f},  // K[2][5]
    {       4.596f,     -10.5183f,      47.3589f,      10.2677f,       12.627f,     -44.8213f},  // K[2][6]
    {    0.379881f,    -0.389245f,      2.07322f,     0.405382f,     0.462157f,    -0.912423f},  // K[2][7]
    {     22.9677f,     -23.9798f,     -42.3335f,      21.3771f,      22.5598f,      21.4827f},  // K[2][8]
    {     2.50105f,     -4.03554f,     -3.27464f,      3.73661f,      3.79208f,     0.197367f},  // K[2][9]
    {    0.681807f,      8.17753f,     -3.98993f,     -10.1914f,     0.947507f,      3.71218f},  // K[3][0]
    {     1.13808f,      10.5199f,     -5.94096f,     -13.6433f,      2.08977f,      5.39406f},  // K[3][1]
    {     10.5617f,     -41.9223f,      10.3498f,      51.5766f,     -22.4015f,     -10.2831f},  // K[3][2]
    {     1.74492f,     -8.30697f,      2.42508f,      9.30878f,     -5.43811f,     -2.25218f},  // K[3][3]
    {       4.596f,      47.3589f,     -10.5183f,     -44.8213f,       12.627f,      10.2677f},  // K[3][4]
    {    0.379881f,      2.07322f,    -0.389245f,    -0.912423f,     0.462157f,     0.405382f},  // K[3][5]
    {     2.70925f,     -7.42684f,      12.1661f,      9.94957f,     -26.4818f,     -9.89115f},  // K[3][6]
    {   0.0435913f,    -0.245292f,     0.888189f,     0.625218f,     -2.04902f,    -0.453017f},  // K[3][7]
    {     22.9677f,     -42.3335f,     -23.9798f,      21.4827f,      22.5598f,      21.3771f},  // K[3][8]
    {     2.50105f,     -3.27464f,     -4.03554f,     0.197367f,      3.79208f,      3.73661f}   // K[3][9]
};
float Offset_Fit_Coefficients[3][6] = {
    {  -0.0153625f,     0.181262f, -2.18624e-15f,     -0.21214f,   3.1123e-15f,  2.61088e-15f},  // theta_l_eq
    {  -0.0153625f, -1.15818e-15f,     0.181262f,  1.86941e-15f,  6.21789e-16f,     -0.21214f},  // theta_r_eq
    {    -0.09442f,  3.96862e-15f,  1.26492e-15f, -5.85297e-15f, -1.93881e-15f, -1.89361e-15f}   // theta_b_eq
};

//定义控制矩阵
float u[10];

//定义拟合增益矩阵
float Fitting_K[4][10];

void Offset_Calc(Compensation_Amount_t *comp,float (*Fit_Coefficients)[6],float L_l,float L_r){
 comp->Gravity_Comp_Theta_l = Fit_Coefficients[0][0] + Fit_Coefficients[0][1]*L_l + Fit_Coefficients[0][2]*L_r+Fit_Coefficients[0][3]*(L_l*L_l)+Fit_Coefficients[0][4]*L_l*L_r+Fit_Coefficients[0][5]*(L_r*L_r);
 comp->Gravity_Comp_Theta_r = Fit_Coefficients[1][0] + Fit_Coefficients[1][1]*L_l + Fit_Coefficients[1][2]*L_r+Fit_Coefficients[1][3]*(L_l*L_l)+Fit_Coefficients[1][4]*L_l*L_r+Fit_Coefficients[1][5]*(L_r*L_r);	
 comp->Gravity_Comp_Theta_b = Fit_Coefficients[2][0] + Fit_Coefficients[2][1]*L_l + Fit_Coefficients[2][2]*L_r+Fit_Coefficients[2][3]*(L_l*L_l)+Fit_Coefficients[2][4]*L_l*L_r+Fit_Coefficients[2][5]*(L_r*L_r);	
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
	u[8] = Max_Output(comp->Gravity_Comp_Theta_b  -  	 body->theta,Theta_B_MAX); u[9] =             0 -   body->d_theta;
	
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
				-u[2]*Fitting_K[0][2]-u[3]*Fitting_K[0][3]
				-u[4]*Fitting_K[0][4]-u[5]*Fitting_K[0][5]
				-u[6]*Fitting_K[0][6]-u[7]*Fitting_K[0][7]
				-u[8]*Fitting_K[0][8]-u[9]*Fitting_K[0][9];

				T[1] = 
				-u[0]*Fitting_K[1][0]-u[1]*Fitting_K[1][1]
				-u[2]*Fitting_K[1][2]-u[3]*Fitting_K[1][3]
				-u[4]*Fitting_K[1][4]-u[5]*Fitting_K[1][5]
				-u[6]*Fitting_K[1][6]-u[7]*Fitting_K[1][7]
				-u[8]*Fitting_K[1][8]-u[9]*Fitting_K[1][9];
				
					//轮电机极性
				T[2] = 
				+u[0]*Fitting_K[2][0]+u[1]*Fitting_K[2][1]
				+u[2]*Fitting_K[2][2]+u[3]*Fitting_K[2][3]
				+u[4]*Fitting_K[2][4]+u[5]*Fitting_K[2][5]
				+u[6]*Fitting_K[2][6]+u[7]*Fitting_K[2][7]
				+u[8]*Fitting_K[2][8]+u[9]*Fitting_K[2][9];
					
				T[3] = 
				+u[0]*Fitting_K[3][0]+u[1]*Fitting_K[3][1]
				+u[2]*Fitting_K[3][2]+u[3]*Fitting_K[3][3]
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
				- u[4]*Fitting_K[0][4]*1.2f- u[5]*Fitting_K[0][5]*1.2f
				-	u[6]*Fitting_K[0][6]*1.2f- u[7]*Fitting_K[0][7]*1.2f
				- u[8]*Fitting_K[0][8]*1.4f- u[9]*Fitting_K[0][9]*3.0f;

				T[1] = 
				- 0												 - u[1]*Fitting_K[1][1]
				- 0												 - u[3]*Fitting_K[1][3]
				- u[4]*Fitting_K[1][4]*1.2f- u[5]*Fitting_K[1][5]*1.2f
				-	u[6]*Fitting_K[1][6]*1.2f- u[7]*Fitting_K[1][7]*1.2f
				- u[8]*Fitting_K[1][8]*1.4f- u[9]*Fitting_K[1][9]*3.0f;
				
					//轮电机极性
				T[2] = 
				+	0												 + u[1]*Fitting_K[2][1]
				+ 0												 + u[3]*Fitting_K[2][3]
				+ u[4]*Fitting_K[2][4]*1.2f+ u[5]*Fitting_K[2][5]*1.2f
				+	u[6]*Fitting_K[2][6]*1.2f+ u[7]*Fitting_K[2][7]*1.2f
				+ u[8]*Fitting_K[2][8]*1.4f+ u[9]*Fitting_K[2][9]*3.0f;

				T[3] = 
				+	0												 + u[1]*Fitting_K[3][1]
				+ 0												 + u[3]*Fitting_K[3][3]
				+ u[4]*Fitting_K[3][4]*1.2f+ u[5]*Fitting_K[3][5]*1.2f
				+	u[6]*Fitting_K[3][6]*1.2f+ u[7]*Fitting_K[3][7]*1.2f
				+ u[8]*Fitting_K[3][8]*1.4f+ u[9]*Fitting_K[3][9]*3.0f;
				
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








