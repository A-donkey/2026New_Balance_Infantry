#include "LQR.h"

//定腿长下增益矩阵,仅用于判断极性
/*		
		{  -2.127449f,   -4.563716f,  -25.968608f,   -4.810233f,    8.134357f,    0.400613f,  -36.897427f,   -1.804595f,   56.905946f,    4.972716f},  // T_r_to_b
    {  -2.127449f,   -4.563716f,   25.968608f,    4.810233f,  -36.897427f,   -1.804595f,    8.134357f,    0.400613f,   56.905946f,    4.972716f},  // T_l_to_b
    {   0.892970f,    1.943427f,   -5.706412f,   -0.880264f,    2.534021f,    0.127543f,   15.015058f,    0.744498f,   12.811744f,    1.560160f},  // T_wr_to_r
    {   0.892970f,    1.943427f,    5.706412f,    0.880264f,   15.015058f,    0.744498f,    2.534021f,    0.127543f,   12.811744f,    1.560160f}   // T_wl_to_l
*/
//变腿长下多项式拟合的全部系数
float P[40][6] = 
{ 		
/*
lqr_Q = diag([100, 1, 4000, 1, 3000, 10, 3000, 10, 20000, 50]);
lqr_R = diag([2, 2, 20, 20]);
*/

		{     -3.7383f,        8.138f,      2.68516f,     -4.93237f,     -12.4751f,      6.56344f},  // K[0][0]
    {    -7.78839f,      18.0633f,      7.04894f,     -9.72945f,      -30.027f,      12.2101f},  // K[0][1]
    {     -14.288f,     -39.2779f,     -48.5543f,      50.0963f,     -21.0255f,       86.753f},  // K[0][2]
    {    -1.49669f,     -9.01818f,     -14.6864f,      11.0438f,      -8.1125f,      23.3559f},  // K[0][3]
    {    -3.49732f,      38.9851f,      36.7524f,     -29.6458f,      30.5608f,     -53.1851f},  // K[0][4]
    {   -0.146043f,      1.48781f,      1.50579f,      0.43919f,      2.82656f,     -2.76824f},  // K[0][5]
    {    -41.2576f,     -10.9605f,      58.9099f,      23.5373f,      -49.011f,     -40.7326f},  // K[0][6]
    {    -2.16689f,    -0.460849f,      3.83996f,      1.27395f,     -3.15076f,     -4.97222f},  // K[0][7]
    {     32.4964f,     -68.0041f,      242.869f,      82.3676f,     -22.1258f,     -291.673f},  // K[0][8]
    {     2.57376f,     -7.50635f,      24.7994f,      9.03656f,     -5.23861f,     -26.6422f},  // K[0][9]
    {     -3.7383f,      2.68516f,        8.138f,      6.56344f,     -12.4751f,     -4.93237f},  // K[1][0]
    {    -7.78839f,      7.04894f,      18.0633f,      12.2101f,      -30.027f,     -9.72945f},  // K[1][1]
    {      14.288f,      48.5543f,      39.2779f,      -86.753f,      21.0255f,     -50.0963f},  // K[1][2]
    {     1.49669f,      14.6864f,      9.01818f,     -23.3559f,       8.1125f,     -11.0438f},  // K[1][3]
    {    -41.2576f,      58.9099f,     -10.9605f,     -40.7326f,      -49.011f,      23.5373f},  // K[1][4]
    {    -2.16689f,      3.83996f,    -0.460849f,     -4.97222f,     -3.15076f,      1.27395f},  // K[1][5]
    {    -3.49732f,      36.7524f,      38.9851f,     -53.1851f,      30.5608f,     -29.6458f},  // K[1][6]
    {   -0.146043f,      1.50579f,      1.48781f,     -2.76824f,      2.82656f,      0.43919f},  // K[1][7]
    {     32.4964f,      242.869f,     -68.0041f,     -291.673f,     -22.1258f,      82.3676f},  // K[1][8]
    {     2.57376f,      24.7994f,     -7.50635f,     -26.6422f,     -5.23861f,      9.03656f},  // K[1][9]
    {    0.480434f,     -3.13164f,      6.10765f,      2.39932f,      1.22735f,     -7.56661f},  // K[2][0]
    {     1.17129f,     -6.81507f,      11.1426f,      4.64686f,      4.53696f,     -15.1828f},  // K[2][1]
    {    -10.4125f,     -9.86861f,      40.4205f,      9.59846f,      20.0064f,     -48.9003f},  // K[2][2]
    {    -1.65347f,     -2.82741f,      7.43522f,      2.77259f,      5.01565f,     -7.81926f},  // K[2][3]
    {     2.20057f,      12.7561f,     -6.01704f,     -10.7882f,     -22.3166f,      7.79564f},  // K[2][4]
    {   0.0454322f,     0.852257f,     -0.22915f,    -0.388597f,     -1.83696f,     0.535155f},  // K[2][5]
    {     7.35839f,     -9.28984f,      52.5796f,      9.29419f,      8.92468f,     -53.5797f},  // K[2][6]
    {    0.383615f,    -0.381498f,      2.29436f,     0.420208f,     0.328632f,     -1.12063f},  // K[2][7]
    {     22.8474f,     -23.2081f,     -43.2164f,       20.836f,      21.9068f,      23.4571f},  // K[2][8]
    {     2.74017f,     -4.25002f,     -3.80617f,      3.81641f,      4.06293f,     0.613183f},  // K[2][9]
    {    0.480434f,      6.10765f,     -3.13164f,     -7.56661f,      1.22735f,      2.39932f},  // K[3][0]
    {     1.17129f,      11.1426f,     -6.81507f,     -15.1828f,      4.53696f,      4.64686f},  // K[3][1]
    {     10.4125f,     -40.4205f,      9.86861f,      48.9003f,     -20.0064f,     -9.59846f},  // K[3][2]
    {     1.65347f,     -7.43522f,      2.82741f,      7.81926f,     -5.01565f,     -2.77259f},  // K[3][3]
    {     7.35839f,      52.5796f,     -9.28984f,     -53.5797f,      8.92468f,      9.29419f},  // K[3][4]
    {    0.383615f,      2.29436f,    -0.381498f,     -1.12063f,     0.328632f,     0.420208f},  // K[3][5]
    {     2.20057f,     -6.01704f,      12.7561f,      7.79564f,     -22.3166f,     -10.7882f},  // K[3][6]
    {   0.0454322f,     -0.22915f,     0.852257f,     0.535155f,     -1.83696f,    -0.388597f},  // K[3][7]
    {     22.8474f,     -43.2164f,     -23.2081f,      23.4571f,      21.9068f,       20.836f},  // K[3][8]
    {     2.74017f,     -3.80617f,     -4.25002f,     0.613183f,      4.06293f,      3.81641f}   // K[3][9]
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








