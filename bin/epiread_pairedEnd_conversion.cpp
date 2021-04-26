#include <iostream>
#include <vector>
#include <sstream>      // std::stringstream
#include <array>
#include <unordered_map>

#include <fstream>
using namespace std;

// Based on N. Loyfer Pattern algorithm

string TAB = "\t";
int MAX_PAT_LEN = 300;
struct CpG_record
{
	
	string next_cpg;
	int index;
	CpG_record(){}

	CpG_record(string _next_cpg, int _index) :
              next_cpg(_next_cpg), index(_index) {}
};
unordered_map<string, CpG_record> dictCpG;
struct SNP_record
{
	char ref;
	char alt;
	string next_snp;
	string sp;
	SNP_record(){}

	SNP_record(char _ref, char _alt,string snp,string _sp) :
              ref(_ref), alt(_alt), next_snp(snp), sp(_sp) {}
};
unordered_map<string, SNP_record> dictSNP;


bool DEBUG = false;
bool SNP_data = true; //file is probably not empty, if file is empty-skip adding SMP data

vector <string> line2tokens(string &line);
//void convert_epiread(string genome_cpg);
void convert_epiread(ofstream& merged_epiread);
//int execute(string cmd, string& output);
string vec2string(vector <string> &vec, string coordinates = string());
//void merge_paired_and_print(vector <string> l1, vector <string> l2, string &genome);
void merge_paired_and_print(vector <string> l1, vector <string> l2,ofstream& merged_epiread);

vector <string> line2tokens(string &line) {
    /** Break string line to words (a vector of string tokens) */
    vector <string> result;
    string cell;
    stringstream lineStream(line);
    while (getline(lineStream, cell, '\t')) {
        result.push_back(cell);
    }
    return result;
}

void convert_epiread(ofstream& merged_epiread) {
    /** parse stdin for epiread paired-end file, sorted by name and order of mate
     * Translate them to single-end like epiread format, and output to a file */

	vector <string> row1, row2, row_from_line;

    bool first_in_pair = true;
    for (string line_str; getline(cin, line_str);) { 
		row_from_line = line2tokens(line_str);
        //  row is first out of a couple of rows
        if (first_in_pair) {
        	//cout <<line_str << "first"<<endl;
			row1 = row_from_line;
			first_in_pair = false;
            continue;
        }
		//the current row has no mates - the QNAME of the lines is not the same
		else if (! first_in_pair && row_from_line[1] != row1[1]) {
			//cout <<line_str << "no mate"<<endl;
			row1 = row_from_line; //the current row is the first of the next copoule
			continue;
		}

		//the current row has no mates - the QNAME of the lines is the same but the chr is not
		else if (row_from_line[0] != row1[0]) {
			//cout <<line_str << "same name not same chr"<<endl;
			row1 = row_from_line; //the current row is the first of the next copoule
			continue;
		}
		//the current row has no mates - the QNAME of the lines is the same, and the chr but the strand is not
		else if (row_from_line[3] != row1[3]) {
			//cout <<line_str << "same name not same strand"<<endl;
			row1 = row_from_line; //the current row is the first of the next copoule
			continue;
		}
		//row is the second mate in the couple
		else {
			//cout <<line_str <<" second row"<<endl;

			row2 = row_from_line;
			first_in_pair = true; //next row will be first row
		}
		
        // process couple of lines. write to stdout
		try {
			  merge_paired_and_print(row1, row2,merged_epiread);
		} catch (std::exception &e) {
			cerr << "Failed! exception:" << endl;
			cerr << e.what() << endl;
		}    
	}
}


string vec2string(vector <string> &vec, string coordinates) {
    /** print a epiread-like vector to stdout, tab separated */
    // vec length must be 8, or 6 if there is no SNP.
	int num_of_column=8;
	if (!SNP_data)
		num_of_column=6;
	string str_vec = vec[0] + coordinates;
	//for (int i=1; i<num_of_column; i++)
	 for (int i=1; i<vec.size(); i++)
		str_vec += TAB + vec[i];
	return str_vec;
}

int CpGFirstIndex(string &locus) {
    /**Get the index of the first CpG according to its locus */
    int start_site = 0;
    auto search = dictCpG.find(locus);
    if (search != dictCpG.end()) {
        start_site = (search->second).index;
    } else {
        // This is an internal error - should never happen.
        throw logic_error("Internal Error. Unknown CpG locus: " + locus);
    }
    return start_site;
    
}

int CpGLastLoci(string &chr,string &pos,int length_cpg) {
	//get the locus of the last CpG according to window and length of string
	if (length_cpg==0) return stoi(pos);
	string locus = chr + TAB + pos;
    auto search = dictCpG.find(locus);
	if (search != dictCpG.end()) {
		for (int i=0;i<length_cpg - 1 ;i++) {
			locus = chr + TAB + (search->second).next_cpg;
			search = dictCpG.find(locus);
			if (search == dictCpG.end()) 
				   throw logic_error("Internal Error. Unknown CpG locus: " + locus);
		}
		return stoi((search->second).next_cpg);

	}

   
    // This is an internal error - should never happen.
    throw logic_error("Internal Error. Unknown CpG locus: " + locus);

}

SNP_record FindSNPRecord(string &locus) {
    /** translate CpG index (in range 1,...,28M~) to dictionary */
    SNP_record variant;
    auto search = dictSNP.find(locus);
    if (search != dictSNP.end()) {
        variant = search->second;
    } else {
        // This is an internal error - should never happen.
        throw logic_error("Internal Error. Unknown SNP locus: " + locus);
    }
    return variant;
}


void initializeDictCpG(string cpg)
{

	int cpg_index=1; 
    vector <string> record, next_record;
	ifstream cpg_file(cpg, ios_base::in);
	string line; 
	string next_cpg;
	
	//get first CpG
	getline(cpg_file, line);
	record = line2tokens(line);	
	
	while (getline(cpg_file, line)) {	
		next_record = line2tokens(line);
		next_cpg = next_record[1];
		dictCpG.insert(make_pair(record[0]+TAB+record[1],CpG_record(next_cpg,cpg_index++)));
		record = next_record;
    }
	//last record-no "next cpg"
	dictCpG.insert(make_pair(record[0]+TAB+record[1],CpG_record("",cpg_index++)));

}


void initializeDictSNP(string snp)
{
	ifstream snp_file(snp, ios::in);
	vector <string> record,next_record;
	string line_str; 
	if ( snp_file.peek() == fstream::traits_type::eof() )
	{
		cerr <<"SNP file is empty"<<endl;
		SNP_data =false;
		return;
	}
	//get first line from snp file
	getline(snp_file, line_str);
	record = line2tokens(line_str);	

	string next_snp;
	char ref,alt;
	while (getline(snp_file, line_str)) {	
		next_record = line2tokens(line_str);	
	//	dictSNP.insert(make_pair(record[0]+TAB+record[1],variant.assign( {record[3],record[4]} ));
		next_snp = next_record[1];
		ref = record[3][0];
		alt = record[4][0];
		//dictSNP.insert(make_pair(record[0]+TAB+record[1],SNP_record(record[3],record[4],next_snp)));
		dictSNP.insert(make_pair(record[0]+TAB+record[1],SNP_record(ref,alt,next_snp,record[7])));
		record = next_record;
    }
	//last record-no "next_snp"
	dictSNP.insert(make_pair(record[0]+TAB+record[1],SNP_record(ref,alt,"",record[7])));

}

string convertChar2srting(char c) 
{
	string s="";
	return s+c;
}



SNP_record checkLocus(string &chr, string &pos, string &strand,char &variant)
{ //validate SNP in specific locus. Apply the rule if needed
 
	string debug_data = "";
	string locus = chr + TAB + pos;
	SNP_record snp = FindSNPRecord(locus);	
		
	if ( strand == "+" && (snp.ref  == 'C' || snp.alt == 'C' ) && (variant == 'C' || variant == 'T' ) ) {
			variant = 'Y';
			return snp;
	} 
	if ( strand == "-" && (snp.ref  == 'G' || snp.alt == 'G' )  && (variant == 'G' || variant == 'A' ) ) {
			variant = 'R';
			return snp;
	}
	
	//if SNP is not in the rule, but matches the ref or alt:
	if (snp.ref == variant ||  snp.alt == variant || variant =='N' || variant == '.')
		 return snp;
	debug_data = "(Ref-" +convertChar2srting(snp.ref)+",Alt-"+convertChar2srting(snp.alt)+",SP-"+snp.sp+")"; 

	throw logic_error("SNP: " + convertChar2srting(variant) + debug_data);  
	

}

string GetFinalSNP(vector <string> &line) 
{  // get the final merged SNP 
	string debug_data;
	string cuttent_snp = line[7];
	string final_snp = "0:";
	int snp_length = line[7].length();
	SNP_record snp_rec = checkLocus(line[0],line[6],line[3],cuttent_snp[0]);
	debug_data = "(Ref-" +convertChar2srting(snp_rec.ref)+",Alt-"+convertChar2srting(snp_rec.alt)+",SP-"+snp_rec.sp+")";
	final_snp += convertChar2srting(cuttent_snp[0]);
	if (DEBUG) {
		final_snp += debug_data;
	}		
	if (snp_length ==1)
		return final_snp; 
	
	string next_pos = snp_rec.next_snp;
	string absolute_pos = line[6]; 	
	for (int i=1; i<snp_length; i++) { 
		SNP_record snp_rec = checkLocus(line[0],next_pos,line[3],cuttent_snp[i]);
		final_snp += ":" + to_string(stoi(next_pos) - stoi(absolute_pos)) + ":" + convertChar2srting(cuttent_snp[i]);
		debug_data = "(Ref-" +convertChar2srting(snp_rec.ref)+",Alt-"+convertChar2srting(snp_rec.alt)+",SP-"+snp_rec.sp+")";
		if (DEBUG) 
			final_snp += debug_data;
		next_pos = snp_rec.next_snp;
	}
	return final_snp;
}


string GetFinalSNPWithoutChecking(vector <string> &line) 
{
	string final_snp = "";
	if (line[6]==".")
		return final_snp;
	string debug_data;
	string cuttent_snp = line[7];
	final_snp = "0:";
	int snp_length = line[7].length();
	string locus = line[0] + TAB + line[6];
	SNP_record snp_rec = FindSNPRecord(locus);
	debug_data = "(Ref-" +convertChar2srting(snp_rec.ref)+",Alt-"+convertChar2srting(snp_rec.alt)+",SP-"+snp_rec.sp+")";
	final_snp += convertChar2srting(cuttent_snp[0]);
	if (DEBUG) {
		final_snp += debug_data;
	}		
	if (snp_length ==1)
		return final_snp; 
	
	string next_pos = snp_rec.next_snp;
	string absolute_pos = line[6]; 	
	for (int i=1; i<snp_length; i++) { 		
		locus = line[0] + TAB + next_pos;
		snp_rec = FindSNPRecord(locus);		
		final_snp += ":" + to_string(stoi(next_pos) - stoi(absolute_pos)) + ":" + convertChar2srting(cuttent_snp[i]);
		debug_data = "(Ref-" +convertChar2srting(snp_rec.ref)+",Alt-"+convertChar2srting(snp_rec.alt)+",SP-"+snp_rec.sp+")";
		if (DEBUG) 
			final_snp += debug_data;
		next_pos = snp_rec.next_snp;
	}
	return final_snp;
}

vector<string> mergeSNP(vector<string> l1, vector<string> l2)
{//change snp to desired format, with relative index for each variant in format 0:var1:12:var2:13:var3

	vector <string> returned_snp;
	
	//if one mate has missing value-make it be the first vector
	if (l2[6] == "." && l1[6] != ".") {
		vector <string> tmp = l1;
		l1 = l2;
		l2 = tmp;
	}
	
	//if both mates have SNP in the same variant-make the one with more variants be the first mate
	if (l1[6] == l2[6] && l1[7].length() < l2[7].length()) {
		vector <string> tmp = l1;
		l1 = l2;
		l2 = tmp;
	}
	
	//if both mates are "-" strand, but same position
	if (l1[6] != "." && l2[6] != "." && l1[4] == l2[4] && stoi(l1[6]) > stoi(l2[6])) {
		vector <string> tmp = l1;
		l1 = l2;
		l2 = tmp;
	}
	
		
	//get SNP-length: number of variants for each line:
	string snp1 = l1[7];
	string snp2 = l2[7];
	int snp1_length = l1[7].length();
	int snp2_length = l2[7].length();
	//if both snp has missing values
	if (l1[6] == "." && l2[6] == ".") returned_snp.assign( {".","."} );
	//get SNP data
	//if one read has missing value (".") than use the other value. 
	else if (l1[6] == "." ) 
			returned_snp.assign( {l2[6],GetFinalSNP(l2)} );
	
	//both mates have values in the SNP column, in the same position. l1 has more variants
	else if (l1[6] == l2[6]) {
		//both mates have the same SNP variants
		if (snp1==snp2) 
				returned_snp.assign( {l2[6],GetFinalSNP(l2)} );
		//same position, different variants in both mates, the first mate has longer SNP variant list
		else 
		{
			//check all common variants between mates are the same. If not- put "N" in the relevant position
			for (int i=0; i<snp2_length; i++)
				if (snp1[i] != snp2[i]) 
					l1[7][i] = 'N';
			//same position, different SNP
			int i=0; 
			while (i<l1[7].length()) {
				if (!(l1[7][i] == 'N'))
					break;
				i++;
			}
			
			//mates disagree for all variants-has only "N" for all variants then put "." (missing values)
			if (i==l1[7].length())
				returned_snp.assign( {".","."} );
			//only some of the mates are missing
			else returned_snp.assign( {l1[6],GetFinalSNP(l1)} );

		}
	}
	//different positions, different SNP for both mates, the first position is always l1
	else {
		//create the final string of variants
		string locus = l1[0] + TAB + l1[6];
		SNP_record snp = FindSNPRecord(locus);
		int i=0;
		//get the first SNP position that overlaps between l1 and l2
		while (i<snp1_length)
		{
			i++;
			//the next SNP position is the first position of l2-SNP
			if (snp.next_snp == l2[6] || i==snp1_length)
				break;
			locus = l1[0] + TAB + snp.next_snp;
			snp = FindSNPRecord(locus);
		}	
		
		//if l1 doesn't overlap l2 at all - check if l1 finishes when l2 starts-append l2 to l1
		if (i==snp1_length) {
			//append "." missing values until snp2 starts, and then append snp2
			while (snp.next_snp !=l2[6])
			{
				l1[7] += ".";
				locus = l1[0] + TAB + snp.next_snp;
				snp = FindSNPRecord(locus);
				
			}	
			l1[7] += snp2;
		}
		
		else { // in case where i<snp1_length, snp1 and snp2 overlap, check overlap positions and append all left l2 to l1
			int j=0;
			for (; i<snp1_length; i++,j++ ) {

				//if they overlap, and variants are identical
				if (snp1[i] == snp2[j]) continue;
				//if they overlap, and variants are differnt than each variant need to be checked to the rule 
				char variant1 = snp1[i];
				char variant2 = snp2[j];
				string next_snp=snp.next_snp;
				snp = checkLocus(l1[0],next_snp,l1[3],variant1);
				snp = checkLocus(l1[0],next_snp,l1[3],variant2);
				//  cerr <<"FFFFFFFFFFFFFFFFFFFFFFFFFFF"<<" " << l1[1] <<endl;
				if (variant1 != variant2)
					l1[7][i] = 'N';
				else 
					l1[7][i] = variant1;
			}
			if (j<snp2_length)  //there is l2 variants left, that doesn't overlap
				for (;j<snp2_length; j++) 
					l1[7] += snp2[j];
		}
		returned_snp.assign( {l1[6],GetFinalSNP(l1)} ); 
	}	
	return returned_snp; 
	
}	
	
string add_coordintes (vector<string> &l1)
{//add first and last coordinates in case where there is SNP data
	int last_snp = 0; // works even if there's no SNP data
	string coordinates = "";
	if (SNP_data && l1[7] != ".") //if there's SNP data 
	{ 	//get last index of SNP
		string str_tmp = l1[7].substr(0,l1[7].rfind(":"));
		str_tmp = str_tmp.substr(str_tmp.rfind(":")+1);
		int last_snp = stoi(l1[6])+stoi(str_tmp);
	}

	//get last index of CpG
	//string window = l1[0] + "\t" + l1[4];
	//int index = CpGFirstIndex(window) + l1[5].length() - 1;
	int last_CpG = CpGLastLoci(l1[0],l1[4], l1[5].length() - 1);

	//insert values to vector and print
	if (SNP_data  && l1[6] != ".")
		//l1.insert(l1.begin()+1, to_string(min(stoi(l1[4]),stoi(l1[6]))));
		coordinates += TAB+ to_string(min(stoi(l1[4]),stoi(l1[6])));
	else // no SNP data, coordinates depends only on first CpG
		//l1.insert(l1.begin()+1, l1[4]);
		coordinates += TAB+ l1[4] ;
	//l1.insert(l1.begin()+2, to_string(max(last_snp,last_CpG)));
	coordinates += TAB+  to_string(max(last_snp,last_CpG));
	return coordinates;
}
 
//void merge_paired_and_print(vector <string> l1, vector <string> l2, string &genome) {
void merge_paired_and_print(vector <string> l1, vector <string> l2, ofstream& merged_epiread) {
    /*Merge two epiread-formated line into one */

	if (!DEBUG) {
		l1[1] = ".";  
		l2[1] = ".";		
	}	
	
	
	bool flag_SNP_identical = (SNP_data) ? l1[6] == l2[6] && l1[7] == l2[7] : true;
	//if l2 doesn't add any information to the read-length, sequence and SNP data:
	if (l1[4] == l2[4] && l1[5] == l2[5] && ( flag_SNP_identical )) { 
		//there is an snp value on the identical lines
		try {
				if ( SNP_data && l1[6] != "." )
					l1[7] = GetFinalSNP(l1);

				//add_coordintes(l1);
				string coordinates = add_coordintes(l1);
				merged_epiread << vec2string(l1,coordinates) <<endl;
			}
		catch (std::exception &e) {
			if (SNP_data) {
				cout << vec2string(l1) << "\t" << GetFinalSNPWithoutChecking(l1) <<endl;   
				cout << vec2string(l2) << "\t" << GetFinalSNPWithoutChecking(l2) <<endl;   
			}
			else {
				cout << vec2string(l1) << endl;
			}   cout << vec2string(l2) << endl;

		}
		return;
	}

	//cerr <<"IMbefore " << vec2string(l1) << " " <<SNP_data<<endl;;

	
	//swap lines if l1 first site comes before l2 first site (minus strand)
    if (stoi(l1[4]) > stoi(l2[4])) {
        vector <string> tmp = l1;
        l1 = l2;
        l2 = tmp;
    }
	

	string pattern1 = l1[5];
	string pattern2 = l2[5];
	int pattern1_len = pattern1.length();
	int pattern2_len = pattern2.length();	
	
	//0-based file
	//int first_cpg_site1 = stoi(l1[4]);
	//int first_cpg_site2 = stoi(l2[4]);
	
	string window1,window2;
	//window1 = l1[0] + "\t" + to_string(first_cpg_site1) + "\t" + to_string(first_cpg_site1+1);
	//window2 = l2[0] + "\t" + to_string(first_cpg_site2) + "\t" + to_string(first_cpg_site2+1);
	window1 = l1[0] + "\t" + l1[4];
	window2 = l2[0] + "\t" + l2[4]; 
	int first_cpg1,first_cpg2;
	try {
		first_cpg1 = CpGFirstIndex(window1);
		first_cpg2 = CpGFirstIndex(window2);
    }
    catch (std::exception &e) {
		cout << vec2string(l1) << endl;   
		cout << vec2string(l2) << endl;
		return;
    }
	
		
	int last_site = max(first_cpg1 + pattern1_len, first_cpg2 + pattern2_len);
	int overall_len = last_site-first_cpg1;

	string merged_pattern;  // output pattern

    if (overall_len > MAX_PAT_LEN) // sanity check: make sure the two reads are not too far apart
    {
       // throw invalid_argument("invalid pairing. merged read is too long ");
    	string output_error = "Problem with:\n" + l1[0] + "\t" + l1[1] + "\t" + l1[2] + "\t" + l1[3] + "\t" + l1[4] + "\n" + l2[0] + "\t" + l2[1] + "\t" + l2[2] + "\t" + l2[3] + "\t" + l2[4]  ; 
    	cerr <<output_error<<endl;
    	return;
    }


    // init merged_pattern with missing values
    for (int i = 0; i < overall_len; i++)
        merged_pattern += ".";

    // set merged_pattern head with pat1
    for (int i = 0; i < pattern1_len; i++)
        merged_pattern[i] = pattern1[i];
  
    // set pattern2 in the adjusted position
    for (int i = 0; i < pattern2_len; i++) {
        int adj_i = i + first_cpg2 - first_cpg1;
        if (merged_pattern[adj_i] == '.') {   // this site was missing from read1
            merged_pattern[adj_i] = pattern2[i];
        } else if ((pattern2[i] != '.') && (merged_pattern[adj_i] != pattern2[i])) {
            // read1 and read2 disagree. treat this case as a missing value for now
            merged_pattern[adj_i] = 'N';
        }
    } 

	vector<string> merged_snp;
	try {
		if (SNP_data) //SNP file is not empty 
		{
			merged_snp = mergeSNP(l1,l2);
			l1[5] = merged_pattern; 
			l1[6] = merged_snp[0];
			l1[7] = merged_snp[1]; 
			
		}
		//add_coordintes(l1);
		merged_epiread << vec2string(l1,add_coordintes(l1)) <<endl;  
		
    }
    catch (std::exception &e) {
			cout << vec2string(l1) << "\t" << GetFinalSNPWithoutChecking(l1) <<endl;   
			cout << vec2string(l2) << "\t" << GetFinalSNPWithoutChecking(l2) <<endl;   
    }
}



int main(int argc, char **argv) {
    clock_t begin = clock();

    try {
		if (!(argc ==4 || argc ==5))
			throw invalid_argument("Usage: epiread_pairedEnd  GENOME_CPG_FILE(0-based, not comressed)  SNP_FILE output_file [DEBUG]");

		string cpg_file = argv[1];
		if (0 == cpg_file.compare (cpg_file.length() - 2, 2, "gz"))
			throw  invalid_argument("GENOME_CPG_FILE must be not compressed");

        if (argc == 5) 
			DEBUG =true;
		ofstream merged_epiread;
		initializeDictCpG(argv[1]);
		initializeDictSNP(argv[2]);
		merged_epiread.open (argv[3]);
		convert_epiread(merged_epiread);
		merged_epiread.close();

		
    }
    catch (std::exception &e) {
        cerr << "Failed! exception:" << endl;
        cerr << e.what() << endl;
        return 1;
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cerr << elapsed_secs << endl;
    return 0;
}
