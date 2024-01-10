#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include <iomanip>
#define PI 3.1415926
using namespace std;

//基本广播星历块
struct EPHEMERISBLOCK
	//每小时一个卫星对应一个基本星历块
{
	//PRN号 
	int PRN;
	double a0, a1, a2;//时间改正数
	//六个轨道参数
	double IODE, Crs, Deltan, M0;// ORBIT - 1
	double Cuc, e, Cus, SqrtA;// ORBIT - 2
	double Toe, Cic, OMEGA, Cis;// ORBIT - 3
	double i0, Crc, omega, OMEGAdot;// ORBIT - 4
	double IDOT, GpsWeekNumber, L2C, L2P;// ORBIT - 5
	double SatAccuracy, SatHealth, TGD, IODC;// ORBIT - 6
};

struct GPSTIME
{
	int weekno;
	double weekSecond;
};

struct Time
{
	int nYear;
	int nMounth;
	int nDay;
	int nHour;
	int nMinute;
	double dSecond;

	Time(int nYear = -1, int nMounth = 0, int nDay = 0, int nHour = 0, int nMinute = 0, double dSecond = 0) :
		nYear(nYear), nMounth(nMounth), nDay(nDay), nHour(nHour), nMinute(nMinute), dSecond(dSecond)
	{
	}
	friend ostream& operator<<(ostream& out, Time t);
	bool operator==(const Time& t); //==运算符重载
};

bool Time::operator==(const Time& t)
{
	return nYear == t.nYear && nMounth == t.nMounth && nDay == t.nDay && nHour == t.nHour && nMinute == t.nMinute &&
		abs(dSecond - t.dSecond) < 1e6;
}

ostream& operator<<(ostream& out, Time t)
{
	out.setf(ios::fixed);
	out << fixed << setprecision(3) << "year:" << t.nYear << "\tmonth:" << t.nMounth << "\tday:" << t.nDay << "\thour:"
		<< t.nHour << "\tminute:" << t.nMinute << "\tsecond:" << t.dSecond << endl;
	out.unsetf(ios::fixed);
	return out;
}

struct Position
{
	Position(int nPRN = -1, double X = 0, double Y = 0, double Z = 0) : PRN(nPRN), X(X), Y(Y), Z(Z)
{
}
friend ostream& operator<<(ostream& out, Position c);

	double X;
	double Y;
	double Z;
	int PRN;
};//卫星位置

ostream& operator<<(ostream& out, Position c)
{
	out.setf(ios::fixed);
	out << "PRN:" << fixed << setprecision(3) << c.PRN << "\tX:" << c.X << "\tY:" << c.Y << "\tZ:" << c.Z << endl;
	out.unsetf(ios::fixed);
	return out;
}

struct Coordinates
{
	Time time;
	//存储当前时间段内的卫星数据
	vector<Position> coordinate;
	friend ostream& operator<<(ostream& out, Coordinates c);
};

int Calendar2GpsTime(int nYear, int nMounth, int nDay, int nHour, int nMinute, double dSecond, double& WeekSecond);
int ReadBrodcastEphemeris(string strEpheNam, int& EphemerisBlockNum, EPHEMERISBLOCK* m_pGpsEphemeris);

string& replace_all(string& src, const string& old_value, const string& new_value) {
	// 每次重新定位起始位置，防止上轮替换后的字符串形成新的old_value
	for (string::size_type pos(0); pos != string::npos; pos += new_value.length()) {
		if ((pos = src.find(old_value, pos)) != string::npos) {
			src.replace(pos, old_value.length(), new_value);
		}
		else break;
	}
	return src;
}


//读广播星历文件，数据存储于上面定义的指针中
//参数：strEpheNam表示广播星历文件的完整路径
//      EphemerisBlockNum 返回读取到的星历块个数
int ReadBrodcastEphemeris(string strEpheNam, int& EphemerisBlockNum, EPHEMERISBLOCK* m_pGpsEphemeris)
{
	int HeadLineNum = 0;
	int WeekNo;
	double WeekSecond;
	//打开文件
	fstream pfEph;
	pfEph.open(strEpheNam, ios::in);
	if (!pfEph.is_open()) return 0;
	//读入头文件
	string strLine;

	while (1)
	{
		getline(pfEph, strLine);
		HeadLineNum++;
		if (strLine.find("END OF HEADER") != string::npos)
			break;
	}
	//计算星历块数
	int AllNum = 0;
	while (!pfEph.eof())
	{
		getline(pfEph, strLine);//这个位置指针没有复位
		AllNum++;
	}
	//临时读入星历块
	EphemerisBlockNum = (AllNum + 1) / 8;
	GPSTIME * pGpsTime = new GPSTIME[EphemerisBlockNum];
	if (!pGpsTime) return 0;
	//将文件指针调整到数据位置
	pfEph.clear();
	pfEph.seekg(0, ios::beg);
	for (int i = 0; i < HeadLineNum; i++)
		getline(pfEph, strLine);
	//定义读取的参数
	int mPrn;//卫星号PRNo 
	int year, month, day, hour, minute;//卫星钟参考时刻
	double   msecond;
	double   a0, a1, a2;//卫星钟飘参数
	double   IODE, Crs, DeltN, M0;//数据星历发布时间，在轨道径向方向上周期改正正弦的振幅
	double   Cuc, e, Cus, sqrtA;//轨道延迹方向上周期改正余弦振幅 、扁心率、轨道延迹方向上周期改正正弦振幅 、长半轴平方根 
	double   Toe, Cic, OMEGA, Cis;//星历参考时刻、轨道倾角周期改正余弦项振幅、参考时刻升交点赤径主项、轨道倾角周期改正正弦项振幅
	double   i0, Crc, omega, OMEGADOT;//参考时间轨道倾角、在轨道径向方向上周期改正余余弦的振幅、近地点角距、升交点赤径在赤道平面中的长期变化
	double   IDOT, L2C, GPSWeek, L2P;////轨道倾角变化率、gps周
	double   AccuracyofSat, HealthofSat, TGD, IODC;//卫星精度、卫星健康、电离层群迟改正数(读取文件可见为0，所以不予读取)


	for (int i = 0; i < EphemerisBlockNum; i++)
	{
		//读取卫星PRN号，星历参考时间
		getline(pfEph, strLine);
		replace_all(strLine, "D", "e");
		sscanf_s(strLine.c_str(), "%d %d %d %d %d %d %lf %lf %lf %lf", &mPrn, &year, &month, &day, &hour, &minute, &msecond, &a0, &a1, &a2);

		year += 2000;

		WeekNo = Calendar2GpsTime(year, month, day, hour, minute, msecond, WeekSecond);
		pGpsTime[i].weekno = WeekNo;
		pGpsTime[i].weekSecond = WeekSecond;

		getline(pfEph, strLine);
		replace_all(strLine, "D", "e");
		sscanf_s(strLine.c_str(), "%lf %lf %lf %lf", &IODE, &Crs, &DeltN, &M0);

		//读 Cuc,e,Cus,sqrtA
		getline(pfEph, strLine);
		replace_all(strLine, "D", "e");
		sscanf_s(strLine.c_str(), "%lf %lf %lf %lf", &Cuc, &e, &Cus, &sqrtA);

		//Toe,Cic,OMEGA,Cis;
		getline(pfEph, strLine);
		replace_all(strLine, "D", "e");
		sscanf_s(strLine.c_str(), "%lf %lf %lf %lf", &Toe, &Cic, &OMEGA, &Cis);

		//i0,Crc,w,OMEGADOT
		getline(pfEph, strLine);
		replace_all(strLine, "D", "e");
		sscanf_s(strLine.c_str(), "%lf %lf %lf %lf", &i0, &Crc, &omega, &OMEGADOT);

		//IDOT,L2Cod,GPSWeek,L2PCod
		getline(pfEph, strLine);
		replace_all(strLine, "D", "e");
		sscanf_s(strLine.c_str(), "%lf %lf %lf %lf", &IDOT, &L2C, &GPSWeek, &L2P);

		//AccuracyofSat,HealthofSat,TGD,IODC
		getline(pfEph, strLine);
		replace_all(strLine, "D", "e");
		sscanf_s(strLine.c_str(), "%lf %lf %lf %lf", &AccuracyofSat, &HealthofSat,&TGD, &IODC);
		//
		getline(pfEph, strLine);//把最后一个读取掉
		//赋值
		m_pGpsEphemeris[i].PRN = mPrn;
		m_pGpsEphemeris[i].a0 = a0;
		m_pGpsEphemeris[i].a1 = a1;
		m_pGpsEphemeris[i].a2 = a2;

		//&IODE, &Crs, &DeltN, &M0
		m_pGpsEphemeris[i].IODE = IODE;
		m_pGpsEphemeris[i].Crs = Crs;
		m_pGpsEphemeris[i].Deltan = DeltN;
		m_pGpsEphemeris[i].M0 = M0;
		//&Cuc, &e, &Cus, &sqrtA
		m_pGpsEphemeris[i].Cuc = Cuc;
		m_pGpsEphemeris[i].e = e;
		m_pGpsEphemeris[i].Cus = Cus;
		m_pGpsEphemeris[i].SqrtA = sqrtA;
		//Toe,Cic,OMEGA,Cis;
		m_pGpsEphemeris[i].Toe = Toe;
		m_pGpsEphemeris[i].Cic = Cic;
		m_pGpsEphemeris[i].OMEGA = OMEGA;
		m_pGpsEphemeris[i].Cis = Cis;
		//i0,Crc,omega,OMEGADOT
		m_pGpsEphemeris[i].i0 = i0;
		m_pGpsEphemeris[i].Crc = Crc;
		m_pGpsEphemeris[i].omega = omega;
		m_pGpsEphemeris[i].OMEGAdot = OMEGADOT;
		//iDOT,L2Cod,GPSWeek,L2PCod
		m_pGpsEphemeris[i].IDOT = IDOT;
		m_pGpsEphemeris[i].L2C = L2C;
		m_pGpsEphemeris[i].L2P = L2P;
		m_pGpsEphemeris[i].GpsWeekNumber = GPSWeek;
		//AccuracyofSat,HealthofSat,TGD,IODC
		m_pGpsEphemeris[i].SatAccuracy = AccuracyofSat;
		m_pGpsEphemeris[i].SatHealth = HealthofSat;
		m_pGpsEphemeris[i].TGD = TGD;
		m_pGpsEphemeris[i].IODC = IODC;
	}
	pfEph.close();
	if (pGpsTime) { delete[]pGpsTime; pGpsTime = NULL; }
	return 1;
}





//从年 月 日 转换为GPS 周秒
int Calendar2GpsTime(int nYear, int nMounth, int nDay, int nHour, int nMinute, double dSecond, double& WeekSecond)
{
	int DayofMonth = 0;
	int DayofYear = 0;
	int weekno = 0;
	int dayofweek;
	int m;
	if (nYear < 1980 || nMounth < 1 || nMounth > 12 || nDay < 1 || nDay > 31)  return -1;//如果时间不合法，函数返回-1
	//计算从1980年到当前的前一年的天数
	for (m = 1980; m < nYear; m++)
	{
		if ((m % 4 == 0 && m % 100 != 0) || (m % 400 == 0))
		{
			DayofYear += 366;
		}
		else
			DayofYear += 365;
	}
	//计算当前一年内从元月到当前前一月的天数
	for (m = 1; m < nMounth; m++)
	{
		if (m == 1 || m == 3 || m == 5 || m == 7 || m == 8 || m == 10 || m == 12)
			DayofMonth += 31;
		else if (m == 4 || m == 6 || m == 9 || m == 11)
			DayofMonth += 30;
		else if (m == 2)
		{
			if ((nYear % 4 == 0 && nYear % 100 != 0) || (nYear % 400 == 0))
				DayofMonth += 29;
			else
				DayofMonth += 28;

		}
	}
	DayofMonth = DayofMonth + nDay - 6;//加上当月天数/减去1980年元月的6日		
	weekno = (DayofYear + DayofMonth) / 7;//计算GPS周
	dayofweek = (DayofYear + DayofMonth) % 7;
	//计算GPS 周秒时间
	WeekSecond = dayofweek * 86400 + nHour * 3600 + nMinute * 60 + dSecond;

	return weekno;
}


//计算卫星位置
Position GetPosition(int year, int month, int day, int hour, int minute, double second, int PRN, int EphemerisBlockNum, EPHEMERISBLOCK* m_pGpsEphemeris) {
	//首先找到卫星在数组中的位置
	double WeekSecond;
	if (Calendar2GpsTime(year, month, day, hour, minute, second, WeekSecond) == -1) cout << "Time Error!" << endl;
	int t = 0;
	for (; t < EphemerisBlockNum; t++) {
		if (m_pGpsEphemeris[t].PRN == PRN) break;
	}
	//if (t == EphemerisBlockNum) {
		//cout << "No Related Satellite!";
		//exit(-1);
	//}
	EPHEMERISBLOCK* m = m_pGpsEphemeris + t;//m为卫星数据指针
	const double GM = 3.986004415E14;//地球引力常数
	const double omegaE = 7.292115147E-5;//地球自转角速度
	double n = sqrt(GM) / pow(m->SqrtA, 3) + m->Deltan;//卫星星历平均角速度
	double M = m->M0 + n * (WeekSecond - m->Toe);//平近点角
	//计算偏近点角
	double E0 = 0;
	double E = M + m->e * sin(E0);
	while (abs(E - E0) >= 1e-12) {
		E0 = E;
		E = M + m->e * sin(E);
	}
	//计算真近点角（注意判断象限）
	double f;
	double  deno = (cos(E) - m->e);//分母
	double  numer = sqrt(1.0 - pow(m->e, 2)) * sin(E);//分子
	f = atan(numer / deno);
	if (deno >= 0 && numer >= 0) {
		f = atan(numer / deno);
	}
	else if (numer >= 0 && deno < 0) {
		f = PI + atan(numer / deno);
	}
	else if (numer < 0 && deno < 0) {
		f = PI + atan(numer / deno);
	}
	else {
		f = 2 * PI + atan(numer / deno);
	}

	double u0 = f + m->omega;//升交距角
	double r0 = m->SqrtA * m->SqrtA * (1 - m->e * cos(E));//卫星向径
	//计算摄动改正项
	double deltu = m->Cuc * cos(2 * u0) + m->Cus * sin(2 * u0);
	double deltr = m->Crc * cos(2 * u0) + m->Crs * sin(2 * u0);
	double delti = m->Cic * cos(2 * u0) + m->Cis * sin(2 * u0);
	//计算改正后的升交距角、卫星向径和轨道倾角
	double u = u0 + deltu;
	double r = r0 + deltr;
	double i = m->i0 + m->IDOT * (WeekSecond - m->Toe) + delti;
	//卫星在轨道坐标系中位置
	double x = r * cos(u);
	double y = r * sin(u);
	//计算升交点经度
	double L = m->OMEGA + m->OMEGAdot * (WeekSecond - m->Toe) - WeekSecond * omegaE;
	//计算卫星在瞬时地球坐标系下的坐标
	Position Final;
	Final.X = x * cos(L) - y * cos(i) * sin(L);
	Final.Y = x * sin(L) + y * cos(i) * cos(L);
	Final.Z = y * sin(i);
	Final.PRN = PRN;
	return Final;
}


Position& Ephemeris2Coordinate(const Time& time, const EPHEMERISBLOCK& GpsEphemeris) //使用const取址以提高效率
{
	//将给定时间转换为GPS周秒
	double WeekSecond;
	int WeekNumber = Calendar2GpsTime(time.nYear, time.nMounth, time.nDay, time.nHour, time.nMinute, time.dSecond,
		WeekSecond);
	//定义一系列相关变量
	double n_0, n;
	const double mu = 3.986004415e14; //第一步变量
	double M, t = WeekSecond; //第二步
	double E0 = 0, E = 0; //第三步
	double f;
	double u1;
	double r1;
	double delta_u, delta_r, delta_i;
	double u, r, i;
	double x, y;
	double L;
	const double w_e = 7.292115147e-5;
	double X, Y, Z;

	//1)计算卫星运行的平均角速度
	n_0 = sqrt(mu) / pow(GpsEphemeris.SqrtA, 3);
	n = n_0 + GpsEphemeris.Deltan;
	//2)计算t 时刻卫星的平近点角
	M = GpsEphemeris.M0 + n * (t - GpsEphemeris.Toe);
	//3)计算偏近点角
	E = 0;
	do
	{
		E0 = E;
		E = M + GpsEphemeris.e * sin(E0);
	} while (abs(E - E0) > 1e-12);
	//4)计算真近点角
	f = atan(sqrt(1 - pow(GpsEphemeris.e, 2)) * sin(E) / (cos(E) - GpsEphemeris.e));
	double former = sqrt(1 - pow(GpsEphemeris.e, 2)) * sin(E), latter = cos(E) - GpsEphemeris.e;
	if (f < 0)f += PI;
	if (former < 0)f += PI;
	//5)计算升交角距（未经改正的）
	u1 = GpsEphemeris.omega + f;
	//6)计算卫星向径（未经改正的）
	r1 = pow(GpsEphemeris.SqrtA, 2) * (1 - GpsEphemeris.e * cos(E));
	//7)计算摄动改正项
	delta_u = GpsEphemeris.Cuc * cos(2 * u1) + GpsEphemeris.Cus * sin(2 * u1);
	delta_r = GpsEphemeris.Crc * cos(2 * u1) + GpsEphemeris.Crs * sin(2 * u1);
	delta_i = GpsEphemeris.Cic * cos(2 * u1) + GpsEphemeris.Cis * sin(2 * u1);
	//8)进行摄动改正
	u = u1 + delta_u;
	r = r1 + delta_r;
	i = GpsEphemeris.i0 + GpsEphemeris.IDOT * (t - GpsEphemeris.Toe) +
		delta_i;
	//9)计算卫星在轨道平面坐标系中的位置
	x = r * cos(u);
	y = r * sin(u);
	//10) 计算升交点经度L
	L = GpsEphemeris.OMEGA + GpsEphemeris.OMEGAdot * (t - GpsEphemeris.Toe) -
		t * w_e;
	//11) 计算卫星在瞬时地球坐标系下的坐标
	X = x * cos(L) - y * cos(i) * sin(L);
	Y = x * sin(L) + y * cos(i) * cos(L);
	Z = y * sin(i);
	Position result = Position(GpsEphemeris.PRN, X, Y, Z);
	return result;
}



int readpreccisonEphemeris(string strfilename,int& EpheremisBlockNum, Coordinates*&m_pGpsEphemeris) {
	ifstream fin(strfilename);
	if (!fin.is_open()) {
		return 0;
	}
	//数出有多少个时段的数据
	EpheremisBlockNum = 0;

	char buf[1024];
	int size = sizeof(buf);
	for (int i = 0; i < 22; i++) {//读取前二十行文件
		fin.getline(buf, size);
	}
	while (fin.getline(buf, size)) {
		if (buf[0] == '*') {
			++EpheremisBlockNum;
		}
	}
	fin.clear();
	fin.seekg(0, ios::beg);
	//申请内存
	m_pGpsEphemeris = new Coordinates[EpheremisBlockNum];
	for (int i = 0; i < 22; i++) {
		fin.getline(buf, size);
	}
	//读文件

	//定义计时器
	int i = -1;

	int PRN;
	double x, y, z;
	int nYear;
	int nMounth;
	int nDay;
	int nHour; 
	int nMinute;
	double dSecond;
	while (fin.getline(buf, size)) {
		if (buf[0] == 'P' && buf[1] == 'G') {
			sscanf_s(buf + 2, "%d %lf %lf %lf", &PRN, &x, &y, &z);
			Position pt;
			x *= 1000;
			y *= 1000;
			z *= 1000;
			pt.X = x;
			pt.Y= y;
			pt.Z = x;
			pt.PRN = PRN;
			m_pGpsEphemeris[i].coordinate.push_back(pt);
		}
		else if (buf[0] == '*') {
			++i;
			sscanf_s(buf + 1, "%d %d %d %d %d %lf", &nYear, &nMounth, &nDay, &nHour, &nMinute, &dSecond);
			m_pGpsEphemeris[i].time.nYear = nYear;
			m_pGpsEphemeris[i].time.nMounth = nMounth;
			m_pGpsEphemeris[i].time.nDay = nDay;
			m_pGpsEphemeris[i].time.nHour = nHour;
			m_pGpsEphemeris[i].time.nMinute = nMinute;
			m_pGpsEphemeris[i].time.dSecond = dSecond;
		}
	}
	fin.close();
	return 1;
}

int Ephemeris2Coordinate(int nYear, int nMounth, int nDay, int nHour, int nMinute, double dSecond, string filepath,
	Position*& output, int& number)
{
	//将给定时间转换为GPS周秒
	double WeekSecond;
	int WeekNumber = Calendar2GpsTime(nYear, nMounth, nDay, nHour, nMinute, dSecond, WeekSecond);

	EPHEMERISBLOCK* m_pGpsEphemeris = nullptr;
	int result = ReadBrodcastEphemeris(filepath, number, m_pGpsEphemeris);
	if (result == 0)
	{
		delete[] m_pGpsEphemeris;
		m_pGpsEphemeris = nullptr;
		return 0;
	}
	output = new Position[number];

	for (int iterator = 0; iterator < number; ++iterator)
	{
		output[iterator] = Ephemeris2Coordinate(Time(nYear, nMounth, nDay, nHour, nMinute, dSecond),m_pGpsEphemeris[iterator]);
	}
	delete[] m_pGpsEphemeris;
	m_pGpsEphemeris = nullptr;
	return 1;
}

ostream& operator<<(ostream& out, Coordinates c)
{
	out.setf(ios::fixed);
	out << fixed << setprecision(3) << c.time;
	for (auto iterator : c.coordinate)
	{
		cout << iterator;
	}
	out.unsetf(ios::fixed);
	return out;
}

int main() {
	int EphemerisBlockNum;
	EPHEMERISBLOCK* m_pGpsEphemeris = NULL;
	m_pGpsEphemeris = new EPHEMERISBLOCK[302];
	ReadBrodcastEphemeris("wdc62540.22n", EphemerisBlockNum, m_pGpsEphemeris);
	cout << "Please input the parameters:(year, month, day, hour, minute, second, PRN)" << endl;
	int year, month, day, hour, minute;
	double second;
	int PRN;
	cin >> year >> month >> day >> hour >> minute >> second >> PRN;
	Position Final = GetPosition(year, month, day, hour, minute, second, PRN, EphemerisBlockNum, m_pGpsEphemeris);
	cout << "The satellite position is:" << endl;
	cout << "X: "; printf_s("%.10f km\n", Final.X / 1000);
	cout << "Y: "; printf_s("%.10f km\n", Final.Y / 1000);
	cout << "Z: "; printf_s("%.10f km\n", Final.Z / 1000);
	delete[]m_pGpsEphemeris;
}

//	EPHEMERISBLOCK* m_pGpsEphemeris1 = new EPHEMERISBLOCK[1];
//	m_pGpsEphemeris1->Cic = 0.130385160446E-07;
//	m_pGpsEphemeris1->Cis = 0.949949026108E-07;
//	m_pGpsEphemeris1->Crc = 0.201875000000E+03;
//	m_pGpsEphemeris1->Crs = 0.406250000000E+01;
//	m_pGpsEphemeris1->Cuc = 0.189989805222E-06;
//	m_pGpsEphemeris1->Cus = 0.912137329578E-05;
//	m_pGpsEphemeris1->Deltan = 0.451411660250E-08;
//	m_pGpsEphemeris1->e = 0.678421219345E-02;
//	m_pGpsEphemeris1->GpsWeekNumber = 931;
//	m_pGpsEphemeris1->i0 = 0.958512160302E-02;
//	m_pGpsEphemeris1->IDOT = -0.253939149013E-09;
//	//double   IDOT, L2C, GPSWeek, L2P;////轨道倾角变化率、？？、gps周
//	//double   AccuracyofSat, HealthofSat, TGD, IODC;//卫星精度、卫星健康、电离层群迟改正数
//
//	m_pGpsEphemeris1->IODE = 0.97E02;
//	m_pGpsEphemeris1->M0 = -0.290282040486E+00;
//	m_pGpsEphemeris1->OMEGA = -0.137835982556E+01;
//	m_pGpsEphemeris1->omega = -0.258419417299E+01;
//	m_pGpsEphemeris1->OMEGAdot = -0.819426989566E-08;
//	m_pGpsEphemeris1->PRN = 1;
//	m_pGpsEphemeris1->SqrtA = 0.515365263176E+04;
//	m_pGpsEphemeris1->Toe = 0.72E+04;
//	Position Final = GetPosition(year, month, day, hour, minute, second, PRN, 1, m_pGpsEphemeris1);
//
//	cout << "The satellite position is:" << endl;
//	cout << "X: " << Final.X << endl;
//	cout << "Y: " << Final.Y << endl;
//	cout << "Z: " << Final.Z << endl;
//}
