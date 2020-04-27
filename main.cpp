#include <cstdio>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <cstring>
#include <thread>

using namespace std;
#define MAX_LEN 7
#define MIN_LEN 3
#define MAX_EDGES   280000  // maximal number of edges 
#define MAX_NODES   2*MAX_EDGES
#define MAX_IN  50  // maximal in degree
#define MAX_OUT 50  // maximal out degree

//data array
char data_buf[MAX_EDGES*3*10];
int32_t ids[MAX_NODES];
int32_t idMap_lg[MAX_NODES];  // map from local id to global id
char strTable[MAX_NODES*10];   // map from local id to glocal id string 

// graph array
int32_t dg[MAX_NODES][MAX_IN+1]; // the first element of each row is for storing the row length
int32_t rdg[MAX_NODES][MAX_OUT+1];
char blocked1[MAX_NODES];   // for thread 1
char blocked2[MAX_NODES];   // for thread 2
char blocked3[MAX_NODES];   // for thread 3
char blocked4[MAX_NODES];   // for thread 4


// result array
char ans3[3*500000*20];
char ans4[4*500000*20];
char ans5[5*1000000*20];
char ans6[6*2000000*20];
char ans7[7*3000000*20];
char *Ans1[] = {ans3, ans4, ans5, ans6, ans7};
char *Ans2[] = {&ans3[3*500000*5], &ans4[4*500000*5], &ans5[5*1000000*5], &ans6[6*2000000*5], &ans7[7*3000000*5]};
char *Ans3[] = {&ans3[3*500000*10], &ans4[4*500000*10], &ans5[5*1000000*10], &ans6[6*2000000*10], &ans7[7*3000000*10]};
char *Ans4[] = {&ans3[3*500000*15], &ans4[4*500000*15], &ans5[5*1000000*15], &ans6[6*2000000*15], &ans7[7*3000000*15]};
char **ans_ptr[] = {Ans1, Ans2, Ans3, Ans4};


class HWCodeCraft
{
    int32_t id;     // thread id
	int32_t stack_ptr = 0;  
    int32_t loop_cnt = 0;
    char *blocked;
    int32_t byte_cnt[5] = {0, 0, 0, 0, 0};  
    int32_t nodeStack[7];  // fixed-size stack

public:
    HWCodeCraft(int32_t _id, char * _blocked): id(_id), blocked(_blocked){}
    void findCircuits(int32_t start_id, int32_t end_id, vector<int32_t> &res)
    {
        unordered_map<int32_t, vector<int32_t>> bag2;   // for storing the second layer nodes
		unordered_set<int32_t> bag3;    // for storing the third layer nodes
		for (int i = start_id; i < end_id; ++i) // search loops containing least id from start_id to end_id
        {
			preSearch(i, bag2, bag3); 
            search(0, i, i, bag2, bag3);    // dfs search
			bag2.clear();
			bag3.clear();
        }
        // write the answer's size
		res[0] = loop_cnt;
		for (int i = 1; i < 6; ++i)
		{
			res[i] = byte_cnt[i-1];
		}
    }

private:

    void preSearch(int32_t s, 
                   unordered_map<int32_t, vector<int32_t>> &bag2, 
                   unordered_set<int32_t> &bag3)
    {
        int32_t j, k, l, len1, len2, len3, tmp1, tmp2;
		len1 = rdg[s][0];
        j = 1;
        while (j <= len1)   // find the index of the first id larger than s
        {
            if (rdg[s][j] > s){break;}
            j++;
        }
        for (; j <= len1; ++j)
        {
            tmp1 = rdg[s][j];
            len2 = rdg[tmp1][0];
            k = 1;
            while (k <= len2)   // find the index of the first id larger than s
            {
                if (rdg[tmp1][k] > s){break;}
                k++;
            }
            for (; k <= len2; ++k)
            {
                tmp2 = rdg[tmp1][k];
                if(bag2.find(tmp2) == bag2.end())   // push the second layer node
                {
                    bag2.emplace(tmp2, vector<int32_t>{tmp1});
                }
                else
                {
                    bag2[tmp2].push_back(tmp1);
                }
                len3 = rdg[tmp2][0];
                l = 1;
                while (l <= len3)   // find the index of the first id larger than s
                {
                    if (rdg[tmp2][l] > s){break;}
                    l++;
                }
                bag3.insert(&rdg[tmp2][l], &rdg[tmp2][len3]+1);     // push the third layer nodes
            }
        }
    }

    void search(int depth, int32_t v, int32_t s,
				unordered_map<int32_t, vector<int32_t>> &bag2,
    			unordered_set<int32_t> &bag3)
    {
        nodeStack[stack_ptr++] = v;
        blocked[v] = 1;
		int32_t i, j, k, nId, len1, len2, tmp1, tmp2;
        char *data_ptr = NULL;
        len1 = dg[v][0];
        i = 1;
        while (i <= len1)   // find the index of the first id larger than s
        {
            if (dg[v][i] > s){break;}
            i++;
        }
        // start DFS search
        for (; i <= len1; ++i)
        {
            nId = dg[v][i];
            // the neighbour is currently not accessible
            if (blocked[nId]){continue;}
			// if the neighbour is in the bag
            auto it = bag2.find(nId);
            if(it != bag2.end())
            {
                nodeStack[stack_ptr++] = nId;
                len2 = it->second.size();
                for (j = 0; j < len2; ++j)
				{
					nId = it->second.at(j);
					if (blocked[nId]){continue;}
					loop_cnt++;
                	nodeStack[stack_ptr++] = nId;
                	// push the answer
					data_ptr = ans_ptr[id][depth];  // get the answer array ptr
                	tmp1 = byte_cnt[depth]; // get the current ptr position
                	for (k = 0; k < depth+2; ++k)
                	{
						tmp2 = 10*nodeStack[k];
						do
						{
							data_ptr[tmp1++] = strTable[tmp2++];
						} while(strTable[tmp2]!=0);
						data_ptr[tmp1++] = ',';
                	}
					tmp2 = 10*nodeStack[k];
					do
					{
						data_ptr[tmp1++] = strTable[tmp2++];
					} while(strTable[tmp2]!=0);
					data_ptr[tmp1++] = '\n';
					byte_cnt[depth] = tmp1;
					stack_ptr--;   
            	}
            	stack_ptr--;    
			}
        	nId = dg[v][i];
        	if (depth < MAX_LEN-4 ||        // if current depth is less than 3 or equal to 3 but able to go back to s
           		(depth == MAX_LEN-4 && bag3.find(nId) != bag3.end()))
        	{
            	search(depth+1, nId, s, bag2, bag3);
        	}	  
        }
        blocked[v] = 0;
        stack_ptr--;
    }
};

// This function must only be used in this file
int32_t getNextPos(const char* buf, int32_t pos, int32_t end)
{ 
    while (buf[pos] != ',' && buf[pos] != '\n' && pos < end)
    {
        pos++;
    }
    return pos+1;
}

// this function read data in the ids array
int32_t readData(string file_name)
{
    int32_t i, j, len, tmp1, tmp2, cur_id;
    // open file
    ifstream fin(file_name, ios::binary|ios::in);
    if (!fin.is_open())
    {
        return -1;
    }
    // read text
    len = fin.seekg(0, ios::end).tellg();
    fin.seekg(0, ios::beg).read(data_buf, static_cast<streamsize>(len));
	fin.close();
    // load data
    i = 0;
    j = 0;
    while (i < len)
    {
        tmp1 = atoi(&data_buf[i]);
        i = getNextPos(data_buf, i, len);
        tmp2 = atoi(&data_buf[i]);
        ids[j++] = tmp1;
        ids[j++] = tmp2;
        i = getNextPos(data_buf, i, len);
        i = getNextPos(data_buf, i, len);
    }
    return j;
}

bool writeResult(string file_name, 
				vector<int32_t> & result_cnt1,
				vector<int32_t> & result_cnt2,
				vector<int32_t> & result_cnt3,
				vector<int32_t> & result_cnt4)

{
    int32_t i, byte_cnt;
    int32_t *ptr = NULL;
    // open file
    ofstream fout(file_name, ios::binary|ios::out);
    if (!fout.is_open())
    {
        return false;
    }
    vector<char> buf(20);
    byte_cnt = sprintf(&buf[0], "%d", result_cnt1[0]+result_cnt2[0]+result_cnt3[0]+result_cnt4[0]);
    buf[byte_cnt++] = '\n';
    fout.write(&buf[0], byte_cnt);
    for (i = 0; i < 5; ++i)
    {
        fout.write(ans_ptr[0][i], result_cnt1[i+1]);
		fout.write(ans_ptr[1][i], result_cnt2[i+1]);
		fout.write(ans_ptr[2][i], result_cnt3[i+1]);
		fout.write(ans_ptr[3][i], result_cnt4[i+1]);
    }
    fout.close();
    return true;
}


int main(int argc, char *argv[])
{
    string file_name = "/data/test_data.txt";
    //string file_name = argv[1];
	string output_file = "/projects/student/result.txt";
    clock_t start, end;
    int32_t i, j, k, len, len1, len2, tmp1, tmp2, tmp3, tmp4;

    //--------------------------read data-----------------------------
    start = clock();
    len = readData(file_name);
    if(len == -1)
    {
        cout << "fail to open file." << endl;
        return -1;
    }
    unordered_map<int32_t, int32_t> idMap_gl;   // map from global id to local id
    j = 0;
    for (i = 0; i < len; i+=2)
    {
        tmp1 = ids[i];
        tmp2 = ids[i+1];
        if (idMap_gl.find(tmp1) == idMap_gl.end())
        {
            idMap_gl.emplace(tmp1, -1);
            idMap_lg[j++] = tmp1;
        }
        if (idMap_gl.find(tmp2) == idMap_gl.end())
        {
            idMap_gl.emplace(tmp2, -1);
            idMap_lg[j++] = tmp2;
        }
    }
    // sort the global id
    sort(idMap_lg, &idMap_lg[j-1]+1);
    for (i = 0; i < j; ++i)
    {
        tmp1 = idMap_lg[i];
        idMap_gl[tmp1] = i;
    }
    end = clock();
    cout << "Time consumption for data reading: " << setprecision(6) << ((double)end - start)/CLOCKS_PER_SEC << "s." << endl;

    //-----------------------create directed graph-----------------------------------
	start = clock();
    // initialize the graph
    for (i = 0; i < j; ++i)
    {
        dg[i][0] = 0;
        rdg[i][0] = 0;
    }
    i = 0;
    do {
        tmp1 = idMap_gl[ids[i++]];
        tmp2 = idMap_gl[ids[i++]];
        tmp3 = dg[tmp1][0]+1;
        dg[tmp1][0] = tmp3;
        tmp4 = rdg[tmp2][0]+1;
        rdg[tmp2][0] = tmp4;
        dg[tmp1][tmp3] = tmp2;
        rdg[tmp2][tmp4] = tmp1;
    } while (i < len);
    len = j;
    for (i = 0; i < len; ++i)   // store the id string
    {
        sprintf(&strTable[10*i], "%d", idMap_lg[i]);
    }
    // sort the edges
    for (i = 0; i < len; ++i)
    {
        tmp1 = dg[i][0];
        sort(dg[i]+1, &dg[i][tmp1]+1);
    }
    for (i = 0; i < len; ++i)
    {
        tmp1 = rdg[i][0];
        sort(rdg[i]+1, &rdg[i][tmp1]+1);
    }
    // initialize the block list
	memset(blocked1, 0, len);
    memset(blocked2, 0, len);
    memset(blocked3, 0, len);
    memset(blocked4, 0, len);
    end = clock();
    cout << "Time consumption for constructing directed graph: " << setprecision(6) << ((double)end - start)/CLOCKS_PER_SEC << "s. " << endl;
    
    //-------------------finding all elementary circuits with required lengths------------------------------------//
    start = clock();
    vector<int32_t> res1(6,0), res2(6,0), res3(6,0), res4(6,0);
    HWCodeCraft solver1(0, blocked1);
	HWCodeCraft solver2(1, blocked2);
	HWCodeCraft solver3(2, blocked3);
	HWCodeCraft solver4(3, blocked4);
	thread t1(&HWCodeCraft::findCircuits, &solver1, 0, 0.13*0.42*len, std::ref(res1));
	thread t2(&HWCodeCraft::findCircuits, &solver2, 0.13*0.42*len, 0.13*len, std::ref(res2));
	thread t3(&HWCodeCraft::findCircuits, &solver3, 0.13*len, 0.13*1.92*len, std::ref(res3));
	thread t4(&HWCodeCraft::findCircuits, &solver4, 0.13*1.92*len, len, std::ref(res4));

	if(t1.joinable()){t1.join();}
	if(t2.joinable()){t2.join();}
	if(t3.joinable()){t3.join();}
	if(t4.joinable()){t4.join();}

	end = clock();
    cout << res1[0]+res2[0]+res3[0]+res4[0] << " circuits found." << endl;
    cout << "Approximate time for finding all the circuits: " << setprecision(6) << ((double)end - start)/CLOCKS_PER_SEC/4 << "s. " << endl;
    
    //------------------print the results--------------------------------------------//
    start = clock();
    if(!writeResult(output_file, res1, res2, res3, res4))
    {
        cout << "fail to create output file." << endl;
        return -1;
    }
    end = clock();
    cout << "Time consumption for outputing result: " << setprecision(6) << ((double)end - start)/CLOCKS_PER_SEC << "s. " << endl;
    
	return 0;
}


