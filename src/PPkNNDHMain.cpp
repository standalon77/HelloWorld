//============================================================================
// Name        : PPkNNDH.cpp
// Author      : PJS
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include <iostream>
#include "PaillierCrypto.h"
#include <time.h>
#include <vector>
#include <thread>
#include <chrono>	// temp

using namespace std;

void Initialization(in_t* tIn, out_t* tOut);
void PrintResult(double* ResTime, int iFoundBit, short sTrd);

//void PPkNNTrd(PaillierCrypto *oCrypto, sync_t* tSync, out_t *tOut, in_t *tIn, unsigned short idx, thp_t* tRecv, th_t* tSend)
void PPkNNTrd(PaillierCrypto *oCrypto, sync_t* tSync, out_t *tOut, in_t *tIn, unsigned short idx, thp_t* tRecv)
{
	//unsigned char* ucSendBuf = new unsigned char[HED_SIZE+2*ENC_SIZE];
	unique_lock<mutex> ulSync(tSync->m, defer_lock);
	pre_t tPre;
	skle_t tSkle;
	unsigned int *pSyncCnt2, *pSyncCnt1;
	int iSrt, iEnd, iBit, iNum;
	bool bFirst = true;

	#ifdef _DEBUG_THREAD
	printf("<<<  %03d - th Thread start  >>>\n", idx);
	#endif

	oCrypto->PreComputation(&tPre, &tSkle, tOut, tIn, idx);


	// ************************   Squared Distance   ************************ //

	#ifdef _DEBUG_MAIN_1
	ulSync.lock();
    (tSync->c1)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", tSync->c1, idx);
	#endif
    if (tSync->c1 < NUM_MAIN_THREAD) {
    	tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD;});
    }
    else {
    	printf("\n************   Squared Distance   ************\n\n");
		tSync->cv.notify_all();		tSync->c2=0;
		tOut->start = time(NULL);						// compute time
    }
	ulSync.unlock();
	#endif

	// compute Squared Distance
	//oCrypto->SquaredDist(tOut, tIn, &tPre, idx, tRecv, tSend, ucSendBuf);
	oCrypto->SquaredDist(tOut, tIn, &tPre, idx, tRecv);

	// print the result of Squared Distance for debugging
	#ifdef _DEBUG_MAIN_1
	ulSync.lock();
    (tSync->c2)++;
    if (tSync->c2 < NUM_MAIN_THREAD) {
    	tSync->cv.wait(ulSync, [&] {return tSync->c2>=NUM_MAIN_THREAD;});
    }
    else {
    	tOut->end = time(NULL);		tOut->dTime[0] = (double) (tOut->end-tOut->start);
    	printf("[DH-%03d] Time : %f (seconds) \n", idx, tOut->dTime[0]);

		oCrypto->DebugOutMain("N\t\t\t", oCrypto->GetPubKey()->n, idx);
		for (int i=0 ; i<DATA_SIZE ; i++) {
			printf("%02d - ", i);
			oCrypto->DebugDecMain("s'_i (plaintext)\t", &(tOut->cS[i]), idx);
		}
		tSync->cv.notify_all();		tSync->c1=0;
    }
	ulSync.unlock();
	#endif


	// ************************   Secure Bit-Decomposition   ************************ //

	#ifdef _DEBUG_MAIN_1
	ulSync.lock();
	(tSync->c1)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", tSync->c1, idx);
	#endif
	if (tSync->c1 < NUM_MAIN_THREAD) {
		tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD;});
	}
	else {
		printf("\n************   Secure Bit-Decomposition   ************\n\n");
		tSync->cv.notify_all();		tSync->c2=0;
		tOut->start = time(NULL);
	}
	ulSync.unlock();
	#endif

	// compute Secure Bit-Decomposition
	//oCrypto->SecBitDecomp(tOut, &tSkle, &tPre, bFirst, idx, tRecv, tSend, ucSendBuf);
	oCrypto->SecBitDecomp(tOut, &tSkle, &tPre, bFirst, idx, tRecv);

	// print the result of Squared Distance for debugging
	#ifdef _DEBUG_MAIN_1
	ulSync.lock();
	(tSync->c2)++;
	if (tSync->c2 < NUM_MAIN_THREAD) {
		tSync->cv.wait(ulSync, [&] {return tSync->c2>=NUM_MAIN_THREAD;});
	}
	else {
		tOut->end = time(NULL);		tOut->dTime[1] = (double) (tOut->end-tOut->start);
		printf("[DH-%03d] Time : %f (seconds) \n", idx, tOut->dTime[1]);

		for (int i=0 ; i<DATA_SIZE ; i++) {
			printf("%02d - ", i);
			oCrypto->DebugDecMain("bit decomposed s'_i (plaintext) ", &(tOut->cS[i]), idx);

		    for (int j=DATA_SQUARE_LENGTH-1 ; j>=0 ; j--) {
				oCrypto->DebugDecBitMain(&(tOut->cSbit[i][j]));
				if (j%8==0)		printf("\t");
				if (j%64==0)	printf("\n");
			}
		}
		tSync->cv.notify_all();		tSync->c1=0;
	}
	ulSync.unlock();
	#endif


	// **********************************   SkLE_s    ********************************** //

	#ifdef _DEBUG_MAIN_1
	ulSync.lock();
	(tSync->c1)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", tSync->c1, idx);
	#endif
	if (tSync->c1 < NUM_MAIN_THREAD) {
		tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD;});
	}
	else {
		printf("\n************   SkLE_s   ************\n\n");
		tSync->cv.notify_all();		tSync->c2=0;
		tOut->start = time(NULL);
	}
	ulSync.unlock();
	#endif

	///////////////////////// SkLE_s start /////////////////////////

	oCrypto->SkLE_s_Init(tOut, &tSkle, &tPre, bFirst, idx);

	oCrypto->SkLE_s_Main(tOut, &tSkle, &tPre, bFirst, idx, tRecv);


//	iNum=-1;
//	for (iBit=DATA_SQUARE_LENGTH ; iBit>=0 ; iBit--) {			// SkLE_s_Main을 위해서 마지막 한번더 실행 필요함.
//		iNum++;
//		#ifdef _DEBUG_THREAD
//		printf("\n<<<  %03d - th bit examination (SkLE_s) >>>\n\n", iBit-1);
//		#endif
//
//		//oCrypto->SkLE_s_Main(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv, tSend, ucSendBuf);
//		oCrypto->SkLE_s_Main(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv);
//
//		ulSync.lock();
//		if (iNum%2 == 0) 		{	pSyncCnt2 = &(tSync->c2);	pSyncCnt1 = &(tSync->c1);	}
//		else 					{	pSyncCnt2 = &(tSync->c1);	pSyncCnt1 = &(tSync->c2);	}
//
//		(*pSyncCnt2)++;
//		if (*pSyncCnt2 < NUM_MAIN_THREAD) {
//			tSync->cv.wait(ulSync, [&] {return *pSyncCnt2>=NUM_MAIN_THREAD;});
//		}
//		else {
//			//oCrypto->SkLE_s_Cmp(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv, tSend, ucSendBuf);
//			oCrypto->SkLE_s_Cmp(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv);
//			tSync->cv.notify_all();		*pSyncCnt1=0;
//		}
//		ulSync.unlock();
//
//		// (tSync->c1), Top-k data를 찾은 경우, 2의 순서를 맞추기 위해서 한번더 수행함.
//		if (tOut->iCmp != 0) {
//			if (iNum%2 == 0) {
//				ulSync.lock();
//				(tSync->c1)++;
//				if (tSync->c1 < NUM_MAIN_THREAD) {
//					tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD;});
//				}
//				else {
//					tSync->cv.notify_all();		tSync->c2=0;
//				}
//				ulSync.unlock();
//			}
//			tSkle.iFoundBit[0] = iBit;
//			break;
//		}
//	}	// (FOR-end) Finding Top-k data 종료
//
//	// 결과 복사
//	iSrt = tOut->iRan[idx];
//	iEnd = tOut->iRan[idx+1];
//	if (tOut->iCmp == 1) {
//		for (int i=iSrt ; i<iEnd ; i++)
//			tOut->cRes[i] = tSkle.cTRes[i-iSrt];
//	}
//	else if (tOut->iCmp == 2) {
//		iNum = tOut->iRan[idx+1]-tOut->iRan[idx];		// iNum 재사용
//		iBit = tOut->iRan[idx];							// iBit 재사용
//		for (int i=0 ; i<iNum ; i++) {
//			paillier_mul(oCrypto->GetPubKey(), &(tOut->cRes[iSrt+i]), &(tSkle.cIRes[i]), &(tSkle.cCan[i]));
//		}
//	}
//	else
//		assert(0);
//
//	///////////////////////// SkLE_s code end /////////////////////////

	// print the result of SkLE_s
	#ifdef _DEBUG_MAIN_1
	ulSync.lock();
    (tSync->c2)++;
    if (tSync->c2 < NUM_MAIN_THREAD) {
    	tSync->cv.wait(ulSync, [&] {return tSync->c2>=NUM_MAIN_THREAD;});
    }
    else {
    	tOut->end = time(NULL);		tOut->dTime[2] = (double) (tOut->end-tOut->start);		// compute time
    	printf("[DH-%03d] Time : %f (seconds) \n", idx, tOut->dTime[2]);

    	printf("[DH-%03d] Found bit : %d \n", idx, tSkle.iFoundBit[0]);
    	for (int i=0 ; i<DATA_SIZE ; i++) {
    		printf("%02d - ", i);
    		oCrypto->DebugDecMain("Res'_i (plaintext) ", &(tOut->cRes[i]), idx);
    	}
		tSync->cv.notify_all();		tSync->c1=0;
    }
	ulSync.unlock();
	#endif


	// ************************   Class Frequency of Top-k Data (CFTKD)   ************************ //
//
//	#ifdef _DEBUG_MAIN_1
//	ulSync.lock();
//	(tSync->c1)++;
//	#ifdef _DEBUG_THREAD
//	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", tSync->c1, idx);
//	#endif
//	if (tSync->c1 < NUM_MAIN_THREAD) {
//		tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD;});
//	}
//	else {
//		printf("\n************   Class Frequency of Top-k Data (CFTKD)   ************\n\n");
//		tSync->cv.notify_all();		tSync->c2=0;
//		tOut->start = time(NULL);						// compute time
//	}
//	ulSync.unlock();
//	#endif
//
//	// compute Class Frequency of Top-k Data (CFTKD)
//	//oCrypto->CFTKD(tOut, tIn, &tPre, idx, tRecv, tSend, ucSendBuf);
//	oCrypto->CFTKD(tOut, tIn, &tPre, idx, tRecv);
//
//	// print the result of CFTKD for debugging
//	#ifdef _DEBUG_MAIN_1
//	ulSync.lock();
//    (tSync->c2)++;
//    if (tSync->c2 < NUM_MAIN_THREAD) {
//    	tSync->cv.wait(ulSync, [&] {return tSync->c2>=NUM_MAIN_THREAD;});
//    }
//    else {
//    	// every thread의 결과 합치기
//		for (int i=0 ; i<CLASS_SIZE ; i++) {
//			mpz_set_ui(tOut->cF[i].c, 1);		// cF[i]=0
//			for (int j=0 ; j<NUM_MAIN_THREAD ; j++)
//				paillier_mul(oCrypto->GetPubKey(), &(tOut->cF[i]), &(tOut->cF[i]), &(tOut->cTF[j][i]));
//		}
//
//		tOut->end = time(NULL);		tOut->dTime[3] = (double) (tOut->end-tOut->start);
//		printf("[DH-%03d] Time : %f (seconds) \n", idx, tOut->dTime[3]);
//
//		// print the result of CFTKD
//		for (int i=0 ; i<CLASS_SIZE ; i++) {
//			printf("%02d - ", i);
//			oCrypto->DebugDecMain("Class Frequency of Top-k Data_i (plaintext)\t", &(tOut->cF[i]), idx);
//		}
//		tSync->cv.notify_all();		tSync->c1=0;
//    }
//	ulSync.unlock();
//	#endif
//
//
//	// ***   두번째 연산을 위한 thread 종료   *** //
//	if (idx >= NUM_MAIN_THREAD2) {
//		printf("(Count: xx) <<<  %03d - th Thread END !!!  >>>\n", idx);
//		return ;
//	}
//	bFirst = false;
//	tSkle.bSkle = false;
//
//	// ************************   Secure Bit-Decomposition for f'_j   ************************ //
//
//	#ifdef _DEBUG_MAIN_1
//	ulSync.lock();
//	(tSync->c1)++;
//	#ifdef _DEBUG_THREAD
//	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", tSync->c1, idx);
//	#endif
//	if (tSync->c1 < NUM_MAIN_THREAD2) {
//		tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD2;});
//	}
//	else {
//		printf("\n************   Secure Bit-Decomposition (f'_j)   ************\n\n");
//		tSync->cv.notify_all();		tSync->c2=0;
//		tOut->start = time(NULL);						// compute time
//	}
//	ulSync.unlock();
//	#endif
//
//	// compute Secure Bit-Decomposition for f'_j
//	//oCrypto->SecBitDecomp(tOut, &tSkle, &tPre, bFirst, idx, tRecv, tSend, ucSendBuf);
//	oCrypto->SecBitDecomp(tOut, &tSkle, &tPre, bFirst, idx, tRecv);
//
//	// print the result of Secure Bit-Decomposition of f'_j
//	#ifdef _DEBUG_MAIN_1
//	ulSync.lock();
//	(tSync->c2)++;
//	if (tSync->c2 < NUM_MAIN_THREAD2) {
//		tSync->cv.wait(ulSync, [&] {return tSync->c2>=NUM_MAIN_THREAD2;});
//	}
//	else {
//		tOut->end = time(NULL);		tOut->dTime[4] = (double) (tOut->end-tOut->start);		// compute time
//		printf("[DH-%03d] Time : %f (seconds) \n", idx, tOut->dTime[4]);
//
//		for (int i=0 ; i<CLASS_SIZE ; i++) {
//			printf("%02d - ", i);
//			oCrypto->DebugDecMain("bit decomposed f'_i (plaintext) ", &(tOut->cF[i]), idx);
//
//		    for (int j=DATA_NUMBER_LENGTH-1 ; j>=0 ; j--) {
//				oCrypto->DebugDecBitMain(&(tOut->cFbit[i][j]));
//				if (j%8==0)		printf("\t");
//				if (j%64==0)	printf("\n");
//			}
//		}
//		tSync->cv.notify_all();		tSync->c1=0;
//	}
//	ulSync.unlock();
//	#endif
//
//
//	// **************************   SkLE_s (k=1, f'_j)    ************************** //
//
//	#ifdef _DEBUG_MAIN_1
//	ulSync.lock();
//	(tSync->c1)++;
//	#ifdef _DEBUG_THREAD
//	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", tSync->c1, idx);
//	#endif
//	if (tSync->c1 < NUM_MAIN_THREAD2) {
//		tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD2;});
//	}
//	else {
//		printf("\n************   IPE-FT1 (f'_j)   ************\n\n");
//		tSync->cv.notify_all();		tSync->c2=0;
//		tOut->start = time(NULL);						// compute time
//	}
//	ulSync.unlock();
//	#endif
//
//	///////////////////////// SkLE_s (k=1, f'_j) code start /////////////////////////
//
//	assert(tOut->iRan2[idx+1] - tOut->iRan2[idx] <= tOut->iRan[idx+1] - tOut->iRan[idx]);
//	assert(NUM_MAIN_THREAD2 <= NUM_MAIN_THREAD);
//
//	oCrypto->SkLE_s_Init(tOut, &tSkle, &tPre, bFirst, idx);
//
//	iNum=-1;
//	for (iBit=DATA_NUMBER_LENGTH ; iBit>=0 ; iBit--) {			// SkLE_s_Main을 위해서 마지막 한번더 실행 필요함.
//		iNum++;
//
//		#ifdef _DEBUG_THREAD
//		printf("\n<<<  %03d - th bit examination (IPE-FT1 for f'_j) >>>\n\n", iBit-1);
//		#endif
//
//		//oCrypto->SkLE_s_Main(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv, tSend, ucSendBuf);
//		oCrypto->SkLE_s_Main(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv);
//
//		ulSync.lock();
//		if (iNum%2 == 0) 	{	pSyncCnt2 = &(tSync->c2);	pSyncCnt1 = &(tSync->c1);	}
//		else 					{	pSyncCnt2 = &(tSync->c1);	pSyncCnt1 = &(tSync->c2);	}
//
//		(*pSyncCnt2)++;
//		if (*pSyncCnt2 < NUM_MAIN_THREAD2) {
//			tSync->cv.wait(ulSync, [&] {return *pSyncCnt2>=NUM_MAIN_THREAD2;});
//		}
//		else {
//			//oCrypto->SkLE_s_Cmp(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv, tSend, ucSendBuf);
//			oCrypto->SkLE_s_Cmp(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv);
//			tSync->cv.notify_all();		*pSyncCnt1=0;
//		}
//		ulSync.unlock();
//
//		// (tSync->c1), Top-k data를 찾은 경우, 2의 순서를 맞추기 위해서 한번더 수행함.
//		if (tOut->iCmp != 0) {
//			if (iNum%2 == 0) {
//				ulSync.lock();
//				(tSync->c1)++;
//				if (tSync->c1 < NUM_MAIN_THREAD2)
//					tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD2;});
//				else {
//					tSync->cv.notify_all();		tSync->c2=0;
//				}
//				ulSync.unlock();
//			}
//			tSkle.iFoundBit[1] = iBit;
//			break;
//		}
//	}	// (FOR-end) Finding Top-k data 종료
//
//	// 결과를 복사
//	iSrt = tOut->iRan2[idx];
//	iEnd = tOut->iRan2[idx+1];
//	if (tOut->iCmp == 1) {
//		for (int i=iSrt ; i<iEnd ; i++)
//			tOut->cFRes[i] = tSkle.cTRes[i-iSrt];
//	}
//	else if (tOut->iCmp == 2) {
//		iNum = tOut->iRan2[idx+1]-tOut->iRan2[idx];		// iNum 재사용
//		iBit = tOut->iRan2[idx];							// iBit 재사용
//		for (int i=0 ; i<iNum ; i++)
//			paillier_mul(oCrypto->GetPubKey(), &(tOut->cFRes[iSrt+i]), &(tSkle.cIRes[i]), &(tSkle.cCan[i]));
//	}
//	else
//		assert(0);
//
//	///////////////////////// SkLE_s (k=1, f'_j) code end /////////////////////////
//
//	// print the result of SkLE_s
//	#ifdef _DEBUG_MAIN_1
//	ulSync.lock();
//    (tSync->c2)++;
//    if (tSync->c2 < NUM_MAIN_THREAD2) {
//    	tSync->cv.wait(ulSync, [&] {return tSync->c2>=NUM_MAIN_THREAD2;});
//    }
//    else {
//    	tOut->end = time(NULL);		tOut->dTime[5] = (double) (tOut->end-tOut->start);		// compute time
//    	printf("[DH-%03d] Time : %f (seconds) \n", idx, tOut->dTime[5]);
//
//    	printf("[DH-%03d] Found bit : %d \n", idx, tSkle.iFoundBit[1]);
//    	for (int i=0 ; i<CLASS_SIZE ; i++) {
//    		printf("%02d - ", i);
//    		oCrypto->DebugDecMain("mf'_i (plaintext) ", &(tOut->cFRes[i]), idx);
//    	}
//		tSync->cv.notify_all();		tSync->c1=0;
//    }
//	ulSync.unlock();
//	#endif
//
//	// checking the allocated size  of GMP
//	#ifdef _DEBUG_Assert
//	assert(tSkle.cTRes[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
//	assert(tSkle.cIRes[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
//	assert(tSkle. cCan[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
////	assert(tSkle.	  cK.c->_mp_alloc == 2*GMP_N_SIZE);
//	#endif
//
//	// ************************   Computing mc'_j   ************************ //
//
//	#ifdef _DEBUG_MAIN_1
//	ulSync.lock();
//    (tSync->c1)++;
//	#ifdef _DEBUG_THREAD
//	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", tSync->c1, idx);
//	#endif
//    if (tSync->c1 < NUM_MAIN_THREAD2) {
//    	tSync->cv.wait(ulSync, [&] {return tSync->c1>=NUM_MAIN_THREAD2;});
//    }
//    else {
//    	printf("\n************   Computing mc'_j   ************\n\n");
//		tSync->cv.notify_all();		tSync->c2=0;
//		tOut->start = time(NULL);						// compute time
//    }
//	ulSync.unlock();
//	#endif
//
//	// compute result class
//	oCrypto->ComputeMc(tOut, tIn, idx);
//
//	// print the result of result class for debugging
//	#ifdef _DEBUG_MAIN_1
//	ulSync.lock();
//    (tSync->c2)++;
//    if (tSync->c2 < NUM_MAIN_THREAD2) {
//    	tSync->cv.wait(ulSync, [&] {return tSync->c2>=NUM_MAIN_THREAD2;});
//    }
//    else {
//    	tOut->end = time(NULL);		tOut->dTime[6] = (double) (tOut->end-tOut->start);		// compute time
//    	printf("[DH-%03d] Time : %f (seconds) \n", idx, tOut->dTime[6]);
//
//		for (int i=0 ; i<CLASS_SIZE ; i++) {
//			//oCrypto->DebugOutMain("mc'_i (ciphertext)", cMc[i]->c, idx);
//			printf("%02d - ", i);
//			oCrypto->DebugDecMain("mc'_i (plaintext)", &(tOut->cMc[i]), idx);
//		}
//		//oCrypto->TerminatePgm(idx, ucSendBuf, tSend, tRecv);
//		short sTrd = oCrypto->TerminatePgm(idx, tRecv);
//		PrintResult(tOut->dTime, tSkle.iFoundBit[0], sTrd);
//		tSync->cv.notify_all();		tSync->c1=0;
//    }
//	ulSync.unlock();
//	#endif

	//delete[] ucSendBuf;

	printf("\n<<<  %03d - th Thread END !!!  >>>\n\n", idx);		// temp

	return;
}

void PrintResult(double* ResTime, int iFoundBit, short sTrd)
{
	// 현재시간 및 시간 결과 출력 (Simple version)
	struct tm *date;
	const time_t t = time(NULL);
	double TotalTime=0;
	int Min, Sec;

	for (int i=0 ; i<7 ; i++)
		TotalTime += ResTime[i];
	date = localtime(&t);

	Min = (int)(TotalTime/60);
	Sec = (int)TotalTime%60;

	printf( "\n[%d/%d/%d %d:%d:%d] %.2f \t %.2f \t %.2f (%d) \t %.2f \t %.2f \t %.2f \t %.2f \t %d min %d sec (%.2f) "
			"\t KeySize=%d \t k=%d \t n=%d \t th=%d,%d  \n" ,
		date->tm_year + 1900 , date->tm_mon + 1 , date->tm_mday , date->tm_hour , date->tm_min , date->tm_sec,
		ResTime[0], ResTime[1], ResTime[2], iFoundBit, ResTime[3], ResTime[4], ResTime[5], ResTime[6], Min, Sec, TotalTime,
		MOD_SIZE, PARAM_K, DATA_SIZE, NUM_MAIN_THREAD, sTrd);

	return ;
}

//void SenderTrd(PaillierCrypto *oCrypto, th_t* tSend)
//{
//	oCrypto->Sender(tSend);
//	return;
//}

void ReceiverTrd(PaillierCrypto *oCrypto, thp_t* tRecv)
{
	oCrypto->Receiver(tRecv);
	return;
}

void Initialization(in_t* tIn, out_t* tOut)
{
	for (int i=0 ; i<CLASS_SIZE ; i++)			tIn->Class[i] = i+1;
	for (int i=0 ; i<DATA_DIMN ; i++)
		#ifdef _DEBUG_INIT_1
		mpz_init(tIn->cQuery[i].c);
		#else
		mpz_init2(tIn->cQuery[i].c, 3*GMP_N_SIZE+1);
		#endif

	for (int i=0 ; i<DATA_SIZE ; i++) {
		for (int j=0 ; j<DATA_DIMN ; j++)
			#ifdef _DEBUG_INIT_1
			mpz_init(tIn->cData[i][j].c);
			#else
			mpz_init2(tIn->cData[i][j].c, 3*GMP_N_SIZE+1);
			#endif
		#ifdef _DEBUG_INIT_1
		mpz_init(tIn->cClass[i].c);
		#else
		mpz_init2(tIn->cClass[i].c, 3*GMP_N_SIZE+1);
		#endif
	}

	for (int i=0 ; i<DATA_SIZE ; i++) {
		for (int j=0 ; j<DATA_SQUARE_LENGTH ; j++)
			#ifdef _DEBUG_INIT_1
			mpz_init(tOut->cSbit[i][j].c);
			#else
			mpz_init2(tOut->cSbit[i][j].c, 2*GMP_N_SIZE*2);
			#endif
		#ifdef _DEBUG_INIT_1
		mpz_inits(tOut->cS[i].c, tOut->cRes[i].c, NULL);
		#else
		mpz_init2(tOut->cS[i].c, 2*GMP_N_SIZE*2);
		mpz_init2(tOut->cRes[i].c, 2*GMP_N_SIZE*2);
		#endif
	}

	for (int i=0 ; i<NUM_MAIN_THREAD ; i++)
		#ifdef _DEBUG_INIT_1
		mpz_init(tOut->cTCnt[i].c);
		#else
		mpz_init2(tOut->cTCnt[i].c, 2*GMP_N_SIZE*2);
		#endif
	#ifdef _DEBUG_INIT_1
			mpz_inits(tOut->cCnt.c, tOut->cCmp.c, NULL);
	#else
	mpz_init2(tOut->cCnt.c, 2*GMP_N_SIZE*2);
	mpz_init2(tOut->cCmp.c, 2*GMP_N_SIZE*2);
	#endif

	for (int i=0 ; i<NUM_MAIN_THREAD ; i++)
		for (int j=0 ; j<CLASS_SIZE ; j++)
			#ifdef _DEBUG_INIT_1
						mpz_init_set_ui(tOut->cTF[i][j].c, 1);		// E_cTF[][] = E(0) : 초기화
			#else
			{ 	mpz_init2(tOut->cTF[i][j].c, 2*GMP_N_SIZE*2);
				mpz_set_ui(tOut->cTF[i][j].c, 1); 			}
			#endif

	// cTRes[j] = f[j] for j=1~v : result of CFTKD
	for (int i=0 ; i<CLASS_SIZE ; i++) {
		for (int j=0 ; j<DATA_NUMBER_LENGTH ; j++)
			#ifdef _DEBUG_INIT_1
			mpz_init(tOut->cFbit[i][j].c);
			#else
			mpz_init2(tOut->cFbit[i][j].c, 2*GMP_N_SIZE*2);
			#endif
		#ifdef _DEBUG_INIT_1
		mpz_inits(tOut->cF[i].c, tOut->cFRes[i].c, tOut->cMc[i].c, NULL);
		#else
		mpz_init2(tOut->cF[i].c, 2*GMP_N_SIZE*2);
		mpz_init2(tOut->cFRes[i].c, 2*GMP_N_SIZE*2);
		mpz_init2(tOut->cMc[i].c, 2*GMP_N_SIZE);
		#endif
	}

	for (int i=0 ; i<7 ; i++)
		tOut->dTime[i] = 0;

	// 각 thread의 데이터 처리 범위 정하기
	#ifdef _DEBUG_Assert
	assert(DATA_SIZE >= NUM_MAIN_THREAD);
	#endif

	int iQuo, iRem;
	iQuo = (int)(DATA_SIZE / NUM_MAIN_THREAD);
	iRem = DATA_SIZE % NUM_MAIN_THREAD;
	tOut->iRan[0] = 0;
	for (int i=0 ; i<NUM_MAIN_THREAD ; i++) {
		tOut->iRan[i+1] = tOut->iRan[i] + iQuo;
		if (i<iRem)
			tOut->iRan[i+1]++;
	}

	#ifdef _DEBUG_Initialization
	std::cout << "[DH] Data range of the 1st Thread : ";
	for (int i=0 ; i<NUM_MAIN_THREAD+1 ; i++)
		std::cout << tOut->iRan[i] << " ";
	printf("\n");
	#endif

	#ifdef _DEBUG_Assert
	assert(tOut->iRan[NUM_MAIN_THREAD]==DATA_SIZE);
	#endif

	// 각 Main thread2의 데이터 처리 범위 정하기
	// NUM_MAIN_THREAD2 = NUM_MAIN_THREAD>=CLASS_SIZE?CLASS_SIZE:NUM_MAIN_THREAD;
	tOut->iRan2[0] = 0;
	if (NUM_MAIN_THREAD >= CLASS_SIZE) {
		for (int i=0 ; i<CLASS_SIZE ; i++)
			tOut->iRan2[i+1] = i+1;
	}
	else {
		iQuo = (int)(CLASS_SIZE / NUM_MAIN_THREAD2);
		iRem = CLASS_SIZE % NUM_MAIN_THREAD2;
		for (int i=0 ; i<NUM_MAIN_THREAD2 ; i++) {
			tOut->iRan2[i+1] = tOut->iRan2[i] + iQuo;
			if (i<iRem)
				tOut->iRan2[i+1]++;
		}
	}

	#ifdef _DEBUG_Initialization
	std::cout << "[DH] Data range of the 2nd Thread : ";
	for (int i=0 ; i<NUM_MAIN_THREAD2+1 ; i++)
		std::cout << tOut->iRan2[i] << " ";
	printf("\n");
	#endif

	#ifdef _DEBUG_Assert
	assert(tOut->iRan2[NUM_MAIN_THREAD2]==CLASS_SIZE);
	#endif

	return;
}

int main() {
//	vector<thread> PPkNNThreads, SenderThreads, ReceiverThreads;
	vector<thread> PPkNNThreads, ReceiverThreads;
	in_t   tIn;
	out_t  tOut;
	sync_t tSync;
//	th_t   tSend;
	thp_t  tRecv;
	PaillierCrypto oCrypto;

	Initialization(&tIn, &tOut);
	tSync.c1 = tSync.c2 = 0;
	for (int i=0 ; i<NUM_MAIN_THREAD ; i++)		tRecv.pa[i] = NULL;

	// distribute public key
	oCrypto.SetPubKey();

	// distribute private key for debugging
	oCrypto.SetPrvKey();

	// input dataset(DO->DH), input query(Q->DH)
	oCrypto.inputDatasetQuery(&tIn);

	// generate PPkNN threads and Communication threads
	for (int i=0 ; i<NUM_MAIN_THREAD ; i++)
		PPkNNThreads.push_back(thread(PPkNNTrd, &oCrypto, &tSync, &tOut, &tIn, (unsigned short)i, &tRecv));
		//PPkNNThreads.push_back(thread(PPkNNTrd, &oCrypto, &tSync, &tOut, &tIn, (unsigned short)i, &tRecv, &tSend));

	//SenderThreads.push_back(thread(SenderTrd, &oCrypto, &tSend));
	ReceiverThreads.push_back(thread(ReceiverTrd, &oCrypto, &tRecv));

	// terminate PPkNN threads and Communication threads

	for (int i=0 ; i<NUM_MAIN_THREAD ; i++)
		PPkNNThreads[i].join();
	std::cout << std::endl << "<<<  (PPkNNDH) Main Thread TERMINATE !!!  >>>"<< std::endl << std::endl;

	//SenderThreads[0].join();
	//std::cout << std::endl << "<<<  (PPkNNDH) Sending Thread TERMINATE !!!  >>>"<< std::endl << std::endl;
	ReceiverThreads[0].join();
	std::cout << std::endl << "<<<  (PPkNNDH) Receiving Thread TERMINATE !!!  >>>"<< std::endl << std::endl;


	// ************************   PPkNN result   ************************ //

	int iTotal=0;
	for (int i=0 ; i<7 ; i++)
		iTotal += tOut.dTime[i];
	std::cout << std::endl << "************   Running Time   ************"<< std::endl << std::endl;
	std::cout << "1. Squared Distance : \t\t\t\t" 				 << tOut.dTime[0] << std::endl;
	std::cout << "2. Secure Bit-Decomposition : \t\t\t" 		 << tOut.dTime[1] << std::endl;
	std::cout << "3. SkLE_s : \t\t\t\t\t" 						 << tOut.dTime[2] << std::endl;
	std::cout << "4. Class Frequency of Top-k Data (CFTKD) : \t" << tOut.dTime[3] << std::endl;
	std::cout << "5. Secure Bit-Decomposition for f'_j : \t\t" 	 << tOut.dTime[4] << std::endl;
	std::cout << "6. IPE-FT1 (f'_j) : \t\t\t\t" 				 << tOut.dTime[5] << std::endl;
	std::cout << "7. Computing mc'_j : \t\t\t\t" 				 << tOut.dTime[6] << std::endl;
	std::cout << "*. Total : \t\t\t\t\t" 				 		 << iTotal		  << std::endl;

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
	assert(tIn.   cQuery[0].c->_mp_alloc == 3*GMP_N_SIZE+1);
	assert(tIn. cData[0][0].c->_mp_alloc == 3*GMP_N_SIZE+1);
	assert(tIn.   cClass[0].c->_mp_alloc == 3*GMP_N_SIZE);
	assert(tOut. 	  cS[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(tOut.cSbit[0][0].c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(tOut.    cRes[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(tOut.   cTCnt[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(tOut.   	   cCnt.c->_mp_alloc <= 2*GMP_N_SIZE*2);
	assert(tOut.	   cCmp.c->_mp_alloc <= 2*GMP_N_SIZE*2);
	assert(tOut.  cTF[0][0].c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(tOut. 	  cF[0].c->_mp_alloc <= 2*GMP_N_SIZE*2);
	assert(tOut.cFbit[0][0].c->_mp_alloc <= 2*GMP_N_SIZE*2);
	assert(tOut.   cFRes[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(tOut. 	 cMc[0].c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	return 0;
}

