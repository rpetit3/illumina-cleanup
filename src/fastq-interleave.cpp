#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <zlib.h>
using namespace std;

int main(int argc,char **argv) {
    ifstream fastq1;
    fastq1.open(argv[1]);
    ifstream fastq2;
    fastq2.open(argv[2]);
    string id1, seq1, plus1, qual1;
    string id2, seq2, plus2, qual2;
    while(true) {
        if(!getline(fastq1,id1,'\n')) break;
        if(!getline(fastq1,seq1,'\n')) break;
        if(!getline(fastq1,plus1,'\n')) break;
        if(!getline(fastq1,qual1,'\n')) break;
        if(!getline(fastq2,id2,'\n')) break;
        if(!getline(fastq2,seq2,'\n')) break;
        if(!getline(fastq2,plus2,'\n')) break;
        if(!getline(fastq2,qual2,'\n')) break;
        cout << id1 << "\n" << seq1 << "\n+\n" << qual1 << endl;
        cout << id2 << "\n" << seq2 << "\n+\n" << qual2 << endl;
    }
    fastq1.close();
    fastq2.close();
    return 0;
}
