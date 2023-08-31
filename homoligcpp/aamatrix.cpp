// a class that handles creation of an amino acid alignment matrix
#include "aamatrix.hpp"

aamatrix::aamatrix(std::string file){
    std::ifstream aafile;
    std::string line,val;
    char delimiter = ' ';
    int k = 0;
    AA_file = file;
    aafile.open(file,std::ifstream::in);
    if (aafile.is_open()){
        while (std::getline(aafile,line)){
            if (line.at(0)=='#') continue;
            std::stringstream sline(line);
            std::vector<float> inner_vec;
            int line_elem=0;
            while (std::getline(sline,val,delimiter)){
                if (val != ""){
                    if (k==0){
                        AA_index[val]=line_elem;
                    } else {
                        if (line_elem>0){
                            inner_vec.push_back(std::stof(val)); // contains every value for the current line in the matrix
                        }
                    }
                    line_elem++;
                }
            }
            if (k>0){
                mat.push_back(inner_vec);
            }
            k++;
        }
    } else {
        std::cout << "file opening error" << std::endl;
    }
    aafile.close();
}

float aamatrix::get_val(std::string AA1, std::string AA2){
    try{
        std::map<std::string,int>::iterator AA1_val = AA_index.find(AA1);
        std::map<std::string,int>::iterator AA2_val = AA_index.find(AA2);
        if (AA1_val != AA_index.end() && AA2_val != AA_index.end()){
            return(mat[AA1_val->second][AA2_val->second]);
        }
    } catch (std::invalid_argument& e){
        std::cerr << e.what() << std::endl;
        std::cerr << AA1 << "or" << AA2 << "are not in " << AA_file << std::endl;
        return(-1);
    }
    return(0);
}

std::vector<std::string> aamatrix::get_kmers(std::string AA_string, int KMER_size){
    int string_length = AA_string.length();
    std::vector<std::string> KMER_vec;
    KMER_vec.reserve(string_length-KMER_size+1);

    for (int i = 0; i+KMER_size <= string_length; i++){
        KMER_vec.push_back(AA_string.substr(i,KMER_size));
    }
    return(KMER_vec);
}


float aamatrix::get_alignment(std::string AA1_string,std::string AA2_string){
    int MIN_SIZE = AA1_string.length();
    int MAX_SIZE = AA2_string.length();
    float MAX_SCORE=0;
    if (AA2_string.length()<AA1_string.length()){
        MIN_SIZE=AA2_string.length();
        MAX_SIZE=AA1_string.length();
    }
    std::vector<std::string> AA1_kmers = get_kmers(AA1_string,MIN_SIZE);
    std::vector<std::string> AA2_kmers = get_kmers(AA2_string,MIN_SIZE);
    for (std::size_t i=0; i<AA1_kmers.size(); i++){
        for (std::size_t j=0; j<AA2_kmers.size(); j++){
            aavector alignment_vector;
            std::string AA1_kmer = AA1_kmers[i];
            std::string AA2_kmer = AA2_kmers[j];
   
            for (std::size_t k=0; k<MIN_SIZE; k++){
            	alignment_vector.add(get_val(std::string(1,AA1_kmer[k]),std::string(1,AA2_kmer[k])));
            }
            float align_val = alignment_vector.getSum()/MAX_SIZE;
            if (align_val>MAX_SCORE){
            	MAX_SCORE=align_val;
                }
        }
        
    }
    
    return(MAX_SCORE);
}
// int aamatrix::check_mat(){
//     std::cout << std::endl;
//     for (int i=0;i<mat.size();i++){
//         for (int j=0;j<mat[i].size();j++){
//             std::cout << mat[i][j] << '\t';
//         }
//         std::cout << std::endl;
//     }
//     std::vector<float> m = mat[0];
//     std::cout << mat[1][1] << std::endl; //why is my amino-acid matrix starting at value (1,1)
//     return(m.size());
// }

// void aamatrix::print_AA_index(void){
//     for (auto it = AA_index.cbegin(); it != AA_index.cend(); ++it){
//         std::cout << it->first << ":" << it->second << "\n";
//     }
// }
