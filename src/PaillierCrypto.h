/*
 * PaillierCrypto.h
 *
 *  Created on: Jun 18, 2019
 *      Author: PJS
 */

#ifndef PAILLIERCRYPTO_H_
#define PAILLIERCRYPTO_H_

#define _DEBUG_INIT_1
#define _DEBUG_MAIN_1
//#define _DEBUG_Assert
#define _DEBUG_THREAD
//#define _DEBUG_Initialization
//#define _DEBUG_SquaredDist
//#define _DEBUG_SecureMultiplication
//#define _DEBUG_SecureBitDecomposition
//#define _DEBUG_EncryptedLSB
//#define _DEBUG_SVR
#define _DEBUG_SkLE_s
#define _DEBUG_Comparison
//#define _DEBUG_CFTKD
//#define _DEBUG_TERMINATE_PGM
//#define _DEBUG_Communication

//불필요한 코드 삭제 및 정리
//send thread 추가

#include <gmp.h>
extern "C"{
	#include "paillier.h"
}
#include "ServerSocket.h"
#include <assert.h>
#include <queue>
#include <mutex>
#include <condition_variable>

const int MOD_SIZE = 1024;						// (check) bit
const int KEY_SIZE = MOD_SIZE/8;				// byte: N
const int ENC_SIZE = KEY_SIZE*2;				// byte: N^2
const int DATA_SQUARE_LENGTH = 12;				// (check) l_1 : data 제곱의 비트수 (실제 데이터 길이)
const int DATA_NUMBER_LENGTH = 5;				// (check) l_2 : data 갯수의 비트수 (실제 데이터 길이)
const int GMP_N_SIZE = MOD_SIZE/64;
const int HED_SIZE = 5;
const int HED_LEN  = 3;

const int DATA_SIZE = 20;						// (check)e
const int DATA_DIMN = 4;							// (check)
const int CLASS_SIZE= 4;							// (check)
const int PARAM_K =   4;							// (check)

const int NUM_MAIN_THREAD  = 1;					// (check) Main thread의 개수.
const int NUM_MAIN_THREAD2 = NUM_MAIN_THREAD>=CLASS_SIZE ? CLASS_SIZE : NUM_MAIN_THREAD;	// Class thread의 개수. (SBD, SkLE_s for f_j)
// 각 thread가 처리하는 데이터의 갯수, (problem) 데이터의 갯수가 정확히 나누어지지 않는 경우 일부 thread는 메모리 낭비가 있을수 있다. (struct로 미리 선언하기 위해 불가피하게 정의함.)
const int NUM_THREAD_DATA  = DATA_SIZE%NUM_MAIN_THREAD==0 ? DATA_SIZE/NUM_MAIN_THREAD : (int)(DATA_SIZE/NUM_MAIN_THREAD)+1;

// Command Tag
const unsigned char COM_MUL1 	= 0x01;
const unsigned char COM_MUL2 	= 0x02;
const unsigned char COM_LSB 	= 0x03;
const unsigned char COM_SVR 	= 0x04;
const unsigned char COM_CMP1 	= 0x05;
const unsigned char COM_CMP2 	= 0x06;
const unsigned char COM_CFTKD 	= 0x07;
const unsigned char COM_TERM 	= 0xFF;

// Main thread의 sync를 관리하기 위한 struct
typedef struct {
	std::mutex m;
	std::condition_variable cv;
	unsigned int c1;
	unsigned int c2;
} sync_t;

//// send를 관리하기 위한 struct
//typedef struct {
//	std::queue<unsigned char*> q;
//	std::mutex m;
//	std::condition_variable cv;
//} th_t;

// receive를 관리하기 위한 struct
typedef struct {
	unsigned char* pa[NUM_MAIN_THREAD];
	std::mutex m;
	std::condition_variable cv;
} thp_t;

typedef struct {
    unsigned int Class[CLASS_SIZE];
    paillier_ciphertext_t cQuery[DATA_DIMN];
    paillier_ciphertext_t cData[DATA_SIZE][DATA_DIMN];		// question: 객체로 된 멤버 변수를 동적생성하는 방법?? (new가 아닌 함수에 의해서 메모리 할당)
    paillier_ciphertext_t cClass[DATA_SIZE];
} in_t;

typedef struct {
	paillier_ciphertext_t cS[DATA_SIZE];						// SSED result
	paillier_ciphertext_t cSbit[DATA_SIZE][DATA_SQUARE_LENGTH];	// SBD result
	paillier_ciphertext_t cRes[DATA_SIZE];						// SkLE_s result
	paillier_ciphertext_t cTCnt[NUM_MAIN_THREAD]; 				// (SkLE_s) thread cnt
	paillier_ciphertext_t cCnt;					 				// (SkLE_s) total cnt
	paillier_ciphertext_t cCmp;					 				// (SkLE_s) comparison
	int					  iCmp;					 				// (SkLE_s) comparison
	paillier_ciphertext_t cTF[NUM_MAIN_THREAD][CLASS_SIZE];		// (CFTKD) thread result
	paillier_ciphertext_t cF[CLASS_SIZE];						// CFTKD result
	paillier_ciphertext_t cFbit[CLASS_SIZE][DATA_NUMBER_LENGTH];// SBD result for f'_j
	paillier_ciphertext_t cFRes[CLASS_SIZE];					// SkLE_s result for f'_j
	paillier_ciphertext_t cMc[CLASS_SIZE];						// computing mc_j
	double dTime[7];
	time_t start;
	time_t end;

    int iRan[NUM_MAIN_THREAD+1];
    int iRan2[NUM_MAIN_THREAD2+1];
} out_t;

typedef struct {
	paillier_plaintext_t  pN1;
	paillier_plaintext_t  pN2;
	paillier_plaintext_t  pL;
	paillier_ciphertext_t c1;
	paillier_ciphertext_t cCls[CLASS_SIZE];
	paillier_ciphertext_t cN1;
	paillier_ciphertext_t cNk;
} pre_t;

typedef struct {
	bool 				  bSkle;									// Cmp func에서 SBD func을 이용하는지 여부를 위한 변수
	paillier_ciphertext_t cCntbit[2][DATA_NUMBER_LENGTH];		// Comparison 함수의 input parameter
	int 				  iKbit		[DATA_NUMBER_LENGTH];		// Comparison 함수의 input parameter
	paillier_ciphertext_t cTRes[NUM_THREAD_DATA];
	paillier_ciphertext_t cIRes[NUM_THREAD_DATA];
	paillier_ciphertext_t cCan[NUM_THREAD_DATA];
	int iFoundBit[2];

	paillier_ciphertext_t cC[NUM_THREAD_DATA];
	paillier_ciphertext_t cK[NUM_THREAD_DATA];
//	unique_lock<mutex> ulSync
} skle_t;

class PaillierCrypto {
private:
	gmp_randstate_t state;
    ServerSocket *mSSocket;
	paillier_prvkey_t* mPrvKey;		// for debugging
    paillier_pubkey_t* mPubKey;

	inline void paillier_ciphertext_from_bytes(paillier_ciphertext_t* ct, void* c, int len);
	inline void paillier_ciphertext_to_bytes(unsigned char* buf, int len, paillier_ciphertext_t* ct );
	inline void SetSendMsg(unsigned char* ucSendPtr, unsigned char* Data, 						 const unsigned short Idx, const unsigned char Tag, const unsigned short Len);
	inline void SetSendMsg(unsigned char* ucSendPtr, unsigned char* Data1, unsigned char* Data2, const unsigned short Idx, const unsigned char Tag, const unsigned short Len1, const unsigned short Len2);
	inline unsigned short Byte2Short(unsigned char* In);
	//inline void Short2Byte(unsigned char* Ret, unsigned short In);
	inline void DebugOut(const mpz_t pData);
	inline void DebugOut(const char* pMsg, 								 const int idx);
	inline void DebugOut(const char* pMsg, const mpz_t pData);
	inline void DebugOut(const char* pMsg, const mpz_t pData, 			 const int idx);
	inline void DebugOut(const char* pMsg, const unsigned char* pData);
	inline void DebugOut(const char* pMsg, short pData, 				 const short idx);
	inline void DebugOut(paillier_pubkey_t* pPubKey);
	inline void DebugDec(paillier_ciphertext_t* cTxt);
	inline void DebugDec(const char* pMsg, paillier_ciphertext_t* cTxt);
	inline void DebugDec(const char* pMsg, paillier_ciphertext_t* cTxt,  const int idx);
	inline void DebugDecBit(paillier_ciphertext_t* cTxt);
	inline void DebugCom(const char* pMsg, char *data, int len, 		 const int idx);
	inline void DebugCom(const char* pMsg, unsigned char *data, int len, const int idx);

//	void SecMul(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cEa, paillier_ciphertext_t* cEb, 	  			   unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	void SecMul(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cEa, 								  			   unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	void EncryptedLSB(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cT, pre_t* tPre, 				  			   unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	int SVR(paillier_ciphertext_t* cEX, paillier_ciphertext_t* cEXi, pre_t* tPre, 						  			   unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	int Comparison(paillier_ciphertext_t* cGamma, paillier_ciphertext_t* cCnt, 							  pre_t *tPre, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);

	void SecMul(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cEa, paillier_ciphertext_t* cEb, 	  			   unsigned short idx, thp_t* tRecv);
	void SecMul(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cEa, 								  			   unsigned short idx, thp_t* tRecv);
	void EncryptedLSB(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cT, 							  pre_t* tPre, unsigned short idx, thp_t* tRecv);
	int SVR(paillier_ciphertext_t* cEX, paillier_ciphertext_t* cEXi, int len, 							  pre_t* tPre, unsigned short idx, thp_t* tRecv);
	int Comparison(paillier_ciphertext_t* cGamma, paillier_ciphertext_t* cCnt,
			paillier_ciphertext_t* cCntbit, paillier_ciphertext_t* cCntCom, int iK, int* iKbit, 		  pre_t *tPre, unsigned short idx, thp_t* tRecv);

public:
	PaillierCrypto();
	virtual ~PaillierCrypto();
	paillier_pubkey_t*			GetPubKey();
	void SetPubKey();
	void SetPrvKey();			// for debugging
	bool inputDatasetQuery(in_t* tIn);
	void PreComputation(pre_t *tPre, skle_t* tSkle, out_t* tOut, in_t* tIn, unsigned short idx);

//	void SquaredDist (out_t* tOut, in_t *tIn, 	pre_t* tPre, 					   unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	void SecBitDecomp(out_t* tOut, 				pre_t* tPre, int iVer,bool bFirst, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	void SkLE_s_Init(out_t* tOut, skle_t* tSkle, pre_t* tPre, 		  bool bFirst, unsigned short idx);
//	void SkLE_s_Main(out_t* tOut, skle_t* tSkle, pre_t* tPre, int bit, bool bFirst, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	void SkLE_s_Cmp (out_t* tOut, skle_t* tSkle, pre_t* tPre, int bit, bool bFirst, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	void CFTKD(out_t* tOut, in_t* tIn, 			pre_t* tPre, 					   unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);
//	void ComputeMc(out_t* tOut, in_t *In, 		   			 					   unsigned short idx);
//	int TerminatePgm(							   			 					   unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr);

	void SquaredDist (out_t* tOut, in_t *tIn,   pre_t* tPre, 					   unsigned short idx, thp_t* tRecv);
	void SecBitDecomp(out_t* tOut, skle_t* tSkle, pre_t* tPre,		  bool bFirst, unsigned short idx, thp_t* tRecv);
	void SkLE_s_Init(out_t* tOut, skle_t* tSkle, pre_t* tPre, 		  bool bFirst, unsigned short idx);
//	void SkLE_s_Main1(out_t* tOut, skle_t* tSkle, pre_t* tPre, int bit, bool bFirst, unsigned short idx, thp_t* tRecv);
	void SkLE_s_Main(out_t* tOut, skle_t* tSkle, pre_t* tPre, bool bFirst, unsigned short idx, thp_t* tRecv);
	void SkLE_s_Cmp (out_t* tOut, skle_t* tSkle, pre_t* tPre, int bit, bool bFirst, unsigned short idx, thp_t* tRecv);
	void CFTKD		 (out_t* tOut, in_t* tIn, 	pre_t* tPre, 					   unsigned short idx, thp_t* tRecv);
	void ComputeMc	 (out_t* tOut, in_t *In, 		   			 				   unsigned short idx);
	int TerminatePgm(							   			 					   unsigned short idx, thp_t* tRecv);

//	int Sender  (th_t*  tSend);
	int Receiver(thp_t* tRecv);
	void DebugOutMain(const char* pMsg, 							 const int idx);
	void DebugOutMain(const char* pMsg, const mpz_t pData, 			 const int idx);
	void DebugDecMain(const char* pMsg, paillier_ciphertext_t* cTxt, const int idx);
	void DebugDecBitMain(paillier_ciphertext_t* cTxt);
};

#endif /* PAILLIERCRYPTO_H_ */
