// COMP9319 2017s1 Assignment 2: Searching BWT Encoded File
// Student Name: Yuxiang Qiu
// Student ID: z5002356
// Implementation:
// 1. For the provided bwt file check its size
// If the size is laege than 2000 bytes then the search will not based on index file
// If not, search will based on index
// 2. For the index based search mode, only need to create index once. If the index file
// size is large than 0, then we know the index has already been created in the previous
// query
// 3. Index file: Since the limitation of memory, we can not open the bwt file and read all
// the content into memory. Instead, we can divide the whole content into k blocks. For each
// block, it records the occurrence from position of start of BWT file of character that is 
// less than ASC II value 128. These information is stored in a int array, the index of array
// corresponds to the ASC value. So assumme there are k blocks, then we got k int array.
// 4. Use backward search to find the FIRST and LAST position in sorted BWT file. Last - First
// is the queries matched.
// 5. For each position between FIRST and LAST do backward and forward concatenation.
// 6. When get a result from previous process it is necessary to check whether the result can
// be matched by second and third query. If all of them are matched, then it is the desired result.
#include<iostream>
#include<fstream>
#include<string>
#include<bitset>
#include<cstring>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <algorithm>

#define Size_Buffer 8196

using namespace std;
void count_c_table_len(int * index, unsigned int &bwt_length){
    for(int i = 0; i < 128; ++i){
        if(index[i] > 0)
            bwt_length += index[i];
    }
}
void count_c_table(int * index, int * c_table, string & str){
    bool count_first_ele_in_c_table = false;
    for(int i = 0;i < 128; ++i){
        if(count_first_ele_in_c_table == false && index[i] > 0){
            c_table[i] = 0;
            str += char(i);
            count_first_ele_in_c_table = true;
        }
        for(int j = 0; j < i; ++j){
            if(index[i] > 0)
                c_table[i] += index[j];  
        }
        if(c_table[i] > 0){
            str += char(i);
        }
    }
}
// get start position of asc_value, the position starts from 1
int get_position(int asc_value, int * index){
    int alphabet_less = 0;
    for(int i = 0; i < 128; ++i){
        if(i == asc_value){
            return alphabet_less;
        }
        else{
            alphabet_less += index[i];
        }
        
    }
    return 0;
}
// get occurrence of asc_value
int occ(int asc_value, int position, FILE * index_file, FILE * bwt_file){
    // need to set the position at start of bwt_file
    fseek(bwt_file, 0, SEEK_SET);
    fseek(index_file, 0, SEEK_SET); 
    int occurrence = 0;
    int block_arr[128];
    char buffer[Size_Buffer];
    memset(buffer,'\0', sizeof(buffer));
    int block_pos = (position + 1) / Size_Buffer;
    // the stop position in buffer
    int buffer_pos = (position + 1) % Size_Buffer;
    // if block_pos is 0 it means position if in the first block so we don not need to locate the block
    if(block_pos > 0){
        // set the position in index file
        fseek(index_file, sizeof(block_arr) * (block_pos - 1), SEEK_SET);
        // read the block before the block
        fread(block_arr,sizeof(block_arr),1,index_file);
        for(int i = 0; i < 128; ++i){
            if(i == asc_value){
                occurrence = block_arr[i];
            }
        }
        // set the position in bwt file
        fseek(bwt_file, sizeof(char) * (block_pos * Size_Buffer), SEEK_SET);
    }
    // read the pos in bwt file
    fread(buffer,sizeof(char),Size_Buffer,bwt_file);
    // convert the asc_value to the corresponding char argument
    char character = (char) asc_value;
    int count = 0;
    for(int i = 0;i < buffer_pos; ++i){
        if(buffer[i] == character){
            ++ count;
        }
    }
    occurrence += count;
    return occurrence;
}
// orientation: ------->
int binary_search(int num_of_blocks, char letter, int index, int * c_table, FILE * index_file, FILE * bwt_file){
    fseek(index_file, 0, SEEK_SET);
    fseek(bwt_file, 0, SEEK_SET);  
    char buffer[Size_Buffer];
    memset(buffer,'\0', sizeof(buffer));
    // need a recorder to store the linear scan in the bwt
    int scan_count = 0;
    // upper and lower bound int array for binary search
    int mid_arr [128] = {};
    // initialize the upper bound with num_of_blocks
    int upper = num_of_blocks;
    // initizlize the mid point
    int mid = num_of_blocks / 2;
    // initialize the lower bound to first block
    int lower = 1;
    int position = 0;
    // first we need to count the order of letter in left hand side(not the bwt)
    int order = index - c_table[int(letter)] + 1;
    int c = 0;
    bool is_previous_diff = false;
    bool is_previous_larger = false;
    int previous_value = 0;
    bool is_previous_same = false;
    int previous_count = 0;
    while(true){
        // read the first block 
        fseek(index_file, sizeof(mid_arr) * (mid - 1), SEEK_SET);
        fread(mid_arr,sizeof(mid_arr),1,index_file);
        // if the info stored in index file is the same with order, it deos not mean we find it. Still need to check whether the previous block
        // is the same with current block
        if(is_previous_same == true){
            if(mid_arr[int(letter)] > order){
                if(mid == 1){
                    fseek(bwt_file, 0, SEEK_SET);
                    fread(buffer,sizeof(char),Size_Buffer,bwt_file);
                }
                else{
                    fseek(index_file, sizeof(mid_arr) * (mid - 2), SEEK_SET);
                    fread(mid_arr,sizeof(mid_arr),1,index_file);
                    position += mid_arr[int(letter)];
                    scan_count += mid_arr[int(letter)];
                    fseek(bwt_file, (mid - 1) * Size_Buffer, SEEK_SET);
                    fread(buffer,sizeof(char),Size_Buffer,bwt_file);
                }
                for(int i = 0;i < Size_Buffer ;++i){
                    ++ position;
                    if(buffer[i] == letter){
                        ++ scan_count; 
                        if(scan_count == order){
                            return position - 1;
                        }
                    }
                }                             
            }
            if(mid_arr[int(letter)] < order){
                if(mid == 1){
                    fseek(bwt_file, Size_Buffer, SEEK_SET);
                    fread(buffer,sizeof(char),Size_Buffer,bwt_file);
                    scan_count += mid_arr[int(letter)];
                    position += Size_Buffer;
                }
                else{
                    fseek(index_file, sizeof(mid_arr) * (mid - 1), SEEK_SET);
                    fread(mid_arr,sizeof(mid_arr),1,index_file);
                    position += Size_Buffer * mid;
                    scan_count = 0;
                    scan_count += mid_arr[int(letter)];
                    fseek(bwt_file, mid * Size_Buffer, SEEK_SET);
                    fread(buffer,sizeof(char),Size_Buffer,bwt_file);
                }
                for(int i = 0;i < Size_Buffer ;++i){
                    ++ position;
                    if(buffer[i] == letter){
                        ++ scan_count; 
                        if(scan_count == order){
                            return position - 1;
                        }
                    }
                }  
            }
            if(mid_arr[int(letter)] == order){

            }
        }
        if(mid_arr[int(letter)] == order){
            // scan the previous blocks till the end to find the first block which has the same order
            while(mid > 1){
                -- mid;
                fseek(index_file, sizeof(mid_arr) * (mid - 1), SEEK_SET);
                fread(mid_arr,sizeof(mid_arr),1,index_file);
                if(mid_arr[int(letter)] == order)
                    continue;
                // if it is not same, then assign the block info to position and scan the next block to get final position
                position += Size_Buffer * mid;
                scan_count += mid_arr[int(letter)];
                fseek(bwt_file, mid * Size_Buffer, SEEK_SET);
                fread(buffer,sizeof(char),Size_Buffer,bwt_file);
                for(int i = 0;i < Size_Buffer ;++i){
                    ++ position;
                    if(buffer[i] == letter){
                        ++ scan_count; 
                        if(scan_count == order){
                            return position - 1;
                        }
                    }
                }                        
            }
            fseek(bwt_file, 0, SEEK_SET);
            fread(buffer,sizeof(char),Size_Buffer,bwt_file);
            scan_count = 0;
            for(int i = 0;i < Size_Buffer ;++i){
                ++ position;
                if(buffer[i] == letter){
                    ++ scan_count; 
                    if(scan_count == order){
                        return position - 1;
                    }
                }
            }                
        }      
        if(mid_arr[int(letter)] > order){
            if(is_previous_same == true)
                continue;
            if(mid == previous_value){
                is_previous_same = true; 
            }
            previous_value = mid;           
            upper = mid;
            mid = (upper + lower) / 2; 
            continue;
        }
        if(mid_arr[int(letter)] < order){
            if(is_previous_same == true)
                continue;
            if(mid == previous_value){
                is_previous_same = true; 
            }
            previous_value = mid; 
            lower = mid;
            mid = (upper + lower) / 2; 
            continue;
        }
        
    }
    return position - 1;
}
string forward_concatenate(string c_table_str, int num_of_blocks, string & query_str, int First, int Last ,FILE * index_file, FILE * bwt_file, int * c_table){
    string concatation;
    char c = query_str[query_str.size() - 1];
    int delta = First - c_table[int(c)] + 1;
    int pos = binary_search(num_of_blocks, query_str.back(), First, c_table, index_file, bwt_file);
    char next_char;
    int j = 0;
    while(true){
        if(next_char == '[')
            break;
        int count = 0;
        for(auto c: c_table_str){
            if(count == c_table_str.size() - 1 && c_table[int(c)] < pos){
                next_char = c;
                concatation += next_char;
                delta = pos - c_table[int(next_char)] + 1;
                pos = c_table[int(c)] + delta - 1;
                pos = binary_search(num_of_blocks, next_char, pos, c_table, index_file, bwt_file);
                break;
            }
            if(c_table[int(c)] == pos){
                next_char = c;
                concatation += next_char;
                delta = pos - c_table[int(next_char)] + 1;
                pos = c_table[int(c)] + delta - 1;
                pos = binary_search(num_of_blocks, next_char, pos, c_table, index_file, bwt_file); 
                break;
            }
            if(c_table[int(c)] > pos){
                next_char = c_table_str[count - 1];
                concatation += next_char;
                delta = pos - c_table[int(next_char)] + 1;
                pos = c_table[int(next_char)] + delta - 1;
                pos = binary_search(num_of_blocks, next_char, pos, c_table, index_file, bwt_file);
                break;
            }
            ++ count;
        }
    }
    concatation.pop_back();
    return concatation;
}
// orientation: <-------
string backward_concatenate(string & query_str, int First, int Last ,FILE * index_file, FILE * bwt_file, int * c_table){
    // need a indicator to check whether the matched result is in square bracket(matched results can not in record_id)
    bool is_in_recordID = true;
    // need a string to store the concatenate info
    string concatation;
    char c = query_str[query_str.size() - 1];
    int previous_pos = c_table[int(c)] + occ(int(c), First, index_file, bwt_file) - 1;   
    int i = 0;
    while(true){
        char buffer[2];
        memset(buffer,'\0', sizeof(buffer));
        fseek(bwt_file, previous_pos, SEEK_SET);
        fread(buffer,sizeof(char),1,bwt_file);
        // concatenate it
        concatation += buffer[0];
        if(buffer[0] == ']'){
            is_in_recordID = false;
        }
        if(buffer[0] == '['){
            break;
        }
        previous_pos = c_table[int(buffer[0])] + occ(int(buffer[0]), previous_pos, index_file, bwt_file) - 1;
        ++ i;
    }
    // reverse the query to do backward search
    reverse(concatation.begin(), concatation.end());
    if(is_in_recordID == true)
        return "";
    else if(is_in_recordID == false)
        return concatation;
    return concatation;
}
int occ_no_idx(int asc_value, int position, FILE * bwt_file){
    // need to set the position at start of bwt_file
    fseek(bwt_file, 0, SEEK_SET);  
    int occurrence = 0;
    char buffer[position + 1];
    memset(buffer,'\0', sizeof(buffer));
    // read the pos in bwt file
    fread(buffer,sizeof(char),position + 1,bwt_file);
    // convert the asc_value to the corresponding char argument
    char character = (char) asc_value;
    int count = 0;
    for(int i = 0;i < position + 1; ++i){
        if(buffer[i] == character){
            ++ count;
        }
    }
    // the block_info + count equals to the position we need
    occurrence += count;  
    return occurrence;
}
int scan_search(char character, int position, string &bwt_content_str, int * c_table, FILE * bwt_file){
    fseek(bwt_file, 0, SEEK_SET); 
    int order = position - c_table[int(character)] + 1;
    int count = 0;
    for(int i = 0;i < bwt_content_str.size(); ++i){
        if(bwt_content_str[i] == character){
            -- order;
            if(order == 0){
                count = i;
                break;
            }      
        }
    }
    return count;
}

string no_idx_backward_concatenate(string & query_str, int First, int Last , FILE * bwt_file, int * c_table){
    // need a indicator to check whether the matched result is in square bracket(matched results can not in record_id)
    bool is_in_recordID = true;
    // need a string to store the concatenate info
    string concatation;
    char c = query_str[query_str.size() - 1];
    int previous_pos = c_table[int(c)] + occ_no_idx(int(c), First, bwt_file) - 1; 
    int i = 0;
    while(true){
        char buffer[2];
        memset(buffer,'\0', sizeof(buffer));
        fseek(bwt_file, previous_pos, SEEK_SET);
        fread(buffer,sizeof(char),1,bwt_file);
        // concatenate it
        concatation += buffer[0];
        if(buffer[0] == ']'){
            is_in_recordID = false;
        }
        if(buffer[0] == '['){
            break;
        }
        previous_pos = c_table[int(buffer[0])] + occ_no_idx(int(buffer[0]), previous_pos, bwt_file) - 1;
        ++ i;
    }
    // reverse the query to do backward search
    reverse(concatation.begin(), concatation.end());
    if(is_in_recordID == true)
        return "";
    else if(is_in_recordID == false)
        return concatation;
    return concatation;   
}

string no_idx_forward_concatenate(string c_table_str, string & query_str, int First, int Last , FILE * bwt_file, int * c_table, string & bwt_content_str){
    string concatation;
    char c = query_str[query_str.size() - 1];
    int delta = First - c_table[int(c)] + 1;
    int pos = scan_search(query_str.back(), First, bwt_content_str, c_table, bwt_file);
    char next_char;
    int j = 0;
    while(true){
        if(next_char == '[')
            break;
        int count = 0;
        for(auto c: c_table_str){
            if(count == c_table_str.size() - 1 && c_table[int(c)] < pos){
                next_char = c;
                concatation += next_char;
                delta = pos - c_table[int(next_char)] + 1;
                pos = c_table[int(c)] + delta - 1;
                pos = scan_search(next_char, pos, bwt_content_str, c_table, bwt_file); 
                break;
            }
            if(c_table[int(c)] == pos){
                next_char = c;
                concatation += next_char;
                delta = pos - c_table[int(next_char)] + 1;
                pos = c_table[int(c)] + delta - 1;
                pos = scan_search(next_char, pos, bwt_content_str, c_table, bwt_file);
                break;
            }
            if(c_table[int(c)] > pos){
                next_char = c_table_str[count - 1];
                concatation += next_char;
                delta = pos - c_table[int(next_char)] + 1;
                pos = c_table[int(next_char)] + delta - 1;
                pos = scan_search(next_char, pos, bwt_content_str, c_table, bwt_file); 
                break;
            }
            ++ count;
        }
    }
    concatation.pop_back();
    return concatation;    
}



int search_without_index(size_t sz, FILE * bwt_file , string &longest_query, string &second_query, string &third_query, int num_of_query){
    string result_str = "";
    unsigned int bwt_length = 0;
    int index[128] = {};
    int c_table[128] = {};
    // store the c_table
    string c_table_str;
    char bwt_content[sz + 1];
    fread(bwt_content,sizeof(char),sz,bwt_file);
    bwt_content[sizeof(bwt_content) - 1] = '\0';
    string bwt_content_str(bwt_content);
    for(auto c: bwt_content_str){
        bitset<8> bits(c);
        int asc_value = bits.to_ulong();
        ++ index[asc_value];
    }
    int First = 0;
    int Last = 0;
    count_c_table(index, c_table, c_table_str);
    // need to count the bwt_length in case the C[current + 1] reach the end of c table
    count_c_table_len(index, bwt_length);
    // reverse the query to do backward search
    reverse(longest_query.begin(), longest_query.end());
    // intialize the backward search algorithm
    int i = 0;
    for(auto c: longest_query){
        if(i == 0){
            First = c_table[int(c)];
            int next_pos;
            for(int loc = 0; loc < c_table_str.size(); ++loc){
                if(loc == c_table_str.size() - 1){
                    next_pos = First + (bwt_length - First) - 1;
                    Last = next_pos;
                    break;
                }
                if(c == c_table_str[loc]){
                    next_pos = int(c_table_str[loc + 1]);
                    Last = c_table[next_pos] - 1;
                    break;
                }
            }
        }
        else{
            First = c_table[int(c)] + occ_no_idx(int(c), First - 1, bwt_file);
            Last = c_table[int(c)] + occ_no_idx(int(c), Last, bwt_file) - 1;
        }
        if(Last < First){
            return 0;
        }
        ++ i;
    }
    string backward_concatenation = "";
    string forward_concatenation = "";
    while(true){
        if(First > Last)
            break;
        int pos = scan_search(longest_query.back(), First, bwt_content_str, c_table, bwt_file);
        backward_concatenation = no_idx_backward_concatenate(longest_query, pos, Last , bwt_file, c_table);
        // check whether backward_concatenation is empty, if it is empty, it means matched result is in square bracket, so it is not we want
        if(backward_concatenation.empty()){
            ++ First;
            result_str = "";
            continue;
        }
        forward_concatenation = no_idx_forward_concatenate(c_table_str, longest_query, First, Last , bwt_file, c_table, bwt_content_str);
        result_str += backward_concatenation;
        result_str += longest_query[longest_query.size() - 1];
        result_str += forward_concatenation;
        string only_text_result = result_str.substr(result_str.find("]") + 1);

        if(num_of_query == 4){
            cout << result_str << endl;  
        }
        if(num_of_query == 5){
            if (only_text_result.find(second_query) != string::npos && second_query != "") {
                cout << result_str << endl;
            }
        }
        if(num_of_query == 6){
            if (only_text_result.find(second_query) != string::npos && 
                only_text_result.find(third_query) != string::npos &&
                second_query != "" &&
                third_query != "") {
                cout << result_str << endl;
            }
        }
        result_str = "";
        ++ First; 
    }
    fclose(bwt_file);
    return 0;
}
int main(int argc, char* argv[]){
    //length of bwt
    unsigned int bwt_length = 0;
    // number of arguments
    int num_of_query = argc;
    // we need to record the amount of block
    int num_of_blocks = 0;
    // whether use index
    bool is_use_idx = true;
    string first_query;
    string second_query;
    string third_query;
    string longest_query;
    if(num_of_query == 4){
        first_query = argv[3];
        longest_query = first_query;
        if (first_query.find("[") != string::npos || first_query.find("]") != string::npos) {
            return 0;
        }
    }
    else if(num_of_query == 5){
        first_query = argv[3];
        second_query = argv[4];
        if (first_query.find("[") != string::npos || 
            first_query.find("]") != string::npos ||
            second_query.find("[") != string::npos ||
            second_query.find("]") != string::npos) {
            return 0;
        }
        if(first_query.size() >= second_query.size())
            longest_query = first_query;
        else{
            longest_query = second_query;
            second_query = first_query;
        }
    }
    else if(num_of_query == 6){
        first_query = argv[3];
        second_query = argv[4];
        third_query =argv[5];
        if (first_query.find("[") != string::npos || 
            first_query.find("]") != string::npos ||
            second_query.find("[") != string::npos ||
            second_query.find("]") != string::npos ||
            third_query.find("[") != string::npos ||
            third_query.find("]") != string::npos) {
            return 0;
        }
        if(first_query.size() >= second_query.size() && first_query.size() >= third_query.size())
            longest_query = first_query;
            string tmp;
            if(second_query.size() <= third_query.size()){
                tmp = second_query;
                second_query = third_query;
                third_query = tmp;
            }
        else if(second_query.size() >= first_query.size() && second_query.size() >= third_query.size()){
            longest_query = second_query;
            if(first_query.size() >= third_query.size()){
                second_query = first_query;
            }
            else{
                second_query = third_query;
                third_query = first_query;
            }
        }
        else if(third_query.size() >= first_query.size() && third_query.size() >= second_query.size()){
            longest_query = third_query;
            string tmp;
            if(first_query.size() >= second_query.size()){
                tmp = second_query;
                second_query = first_query;
                third_query = tmp;
            }
            else{
                third_query = first_query;
            }
        }
    }
    string result_str = "";
    int file_length = 0;
    // used for count for C table
	int index [128] = {};
    // this int array is used for store the block info read from index file
    int block_index[128];
    // C table
    int c_table[128] = {};
	char buffer[Size_Buffer];
    memset(buffer,'\0', sizeof(buffer));
	FILE * bwt_file = fopen(argv[1],"r");
	FILE * index_file = fopen(argv[2],"ab+");

    // check the file size to decide whether use index or not to do the searching
    fseek(bwt_file, 0L, SEEK_END);
    size_t bwt_sz = ftell(bwt_file);
    fseek(bwt_file, 0, SEEK_SET);

    //check the index file size, if it is not empty then no need to create index file again
    //simply jump to the search part
    fseek(index_file, 0L, SEEK_END);
    size_t idx_sz = ftell(index_file);
    fseek(index_file, 0, SEEK_SET);

    if(bwt_sz < 2048){
        int judge = search_without_index(bwt_sz, bwt_file , longest_query, second_query, third_query, num_of_query);
        fclose(index_file);
        return 0;
    }


	while(!feof(bwt_file)){
		fread(buffer,sizeof(char),Size_Buffer,bwt_file);
		for(auto c: buffer){
			if(c == '\0'){
                // only write to index file when it has not created yet
                if(idx_sz == 0){
                    fwrite(index,sizeof(index),1,index_file);
                }   
                ++ num_of_blocks;
                break;
            }
            // if not stop the file length increses by 1
            ++ file_length;
            bitset<8> bits(c);
            int asc_value = bits.to_ulong();
            ++ index[asc_value];
		}
        // only write to index file when it has not created yet
        if(idx_sz == 0){
            fwrite(index,sizeof(index),1,index_file);
        } 
        ++ num_of_blocks;
        memset(buffer,'\0', sizeof(buffer));
	}
    fclose(index_file);

    index_file = NULL;
    index_file = fopen(argv[2],"rb");
    // index file is created start search
    int First = 0;
    int Last = 0;
    // store the c_table
    string c_table_str;
    //count the C table
    count_c_table(index, c_table, c_table_str);
    // need to count the bwt_length in case the C[current + 1] reach the end of c table
    count_c_table_len(index, bwt_length);
    // reverse the query to do backward search
    reverse(longest_query.begin(), longest_query.end());
    // intialize the backward search algorithm
    int i = 0;
    for(auto c: longest_query){
        if(i == 0){
            First = c_table[int(c)];
            int next_pos;
            for(int loc = 0; loc < c_table_str.size(); ++loc){
                if(loc == c_table_str.size() - 1){
                    next_pos = First + (bwt_length - First) - 1;
                    Last = next_pos;
                    break;
                }
                if(c == c_table_str[loc]){
                    next_pos = int(c_table_str[loc + 1]);
                    Last = c_table[next_pos] - 1;
                    break;
                }
            }
        }
        else{
            First = c_table[int(c)] + occ(int(c), First - 1, index_file, bwt_file);
            Last = c_table[int(c)] + occ(int(c), Last, index_file, bwt_file) - 1;
        }
        if(Last < First){
            return 0;
        }
        ++ i;
    }
    string backward_concatenation = "";
    string forward_concatenation = "";
    while(true){
        if(First > Last)
            break;
        int pos = binary_search(num_of_blocks, longest_query.back(), First, c_table, index_file, bwt_file);
        // First is the index we get from the first character in query
        backward_concatenation = backward_concatenate(longest_query, pos, Last ,index_file, bwt_file, c_table);
        // check whether backward_concatenation is empty, if it is empty, it means matched result is in square bracket, so it is not we want
        if(backward_concatenation.empty()){
            ++ First;
            result_str = "";
            continue;
        }
        forward_concatenation = forward_concatenate(c_table_str, num_of_blocks, longest_query, First, Last ,index_file, bwt_file, c_table);
        result_str += backward_concatenation;
        result_str += longest_query[longest_query.size() - 1];
        result_str += forward_concatenation;

        string only_text_result = result_str.substr(result_str.find("]") + 1);

        if(num_of_query == 4){
            cout << result_str << endl;  
        }
        if(num_of_query == 5){
            if (only_text_result.find(second_query) != string::npos && second_query != "") {
                cout << result_str << endl;
            }
        }
        if(num_of_query == 6){
            if (only_text_result.find(second_query) != string::npos && 
                only_text_result.find(third_query) != string::npos &&
                second_query != "" &&
                third_query != "") {
                cout << result_str << endl;
            }
        }
        result_str = "";
        ++ First; 
    } 
	fclose(bwt_file);
	fclose(index_file);
	return 0;
}