#include <iostream>
#include <cassert>

using namespace std;

bool isValidBase(char base){

    bool isValid = false;

    if (base == 'A' || base == 'C' || base == 'G' || base == 'T'){
        isValid = true;
    }
    else{
        isValid = false;
    }
    return isValid;
}

bool isValidStrand(string strand){

    bool isValid = false;
    int length = strand.length();

    for(int i = 0; i < length; i++){
        if(!isValidBase(strand[i])){
            isValid = false;
            break;
        }
        else{
            isValid = true;
        }
    }

    return isValid;
}

double strandSimilarity(string strand1, string strand2){
    double sim = 0;
    double count = 0;

    double length1 = strand1.length();
    double length2 = strand2.length();

    if(length1 != length2){
        return sim;
    }
    else{
        for (int i = 0; i < length1; i++){
            if(strand1[i] == strand2[i]){
                count++;
            }
        }

        sim = count/length1;
        return sim;
    }
    return sim;
}

int bestStrandMatch(string input_strand, string target_strand){
    int startIndex = 0;
    int inputl = input_strand.length();
    int targetl = target_strand.length();
    double finalsim = 0.0;
    double sim = 0.0;

    if(input_strand.length() < target_strand.length()){
        cout << "Best similarity score: 0.0" << endl;
        return -1;
    }
    else{
        for (int i = 0; i < inputl; i++){
            sim = strandSimilarity(input_strand.substr(i, targetl), target_strand);
            if (finalsim < sim){
                finalsim = sim;
                startIndex = i;
            }
        }
    }

    cout << "Best similarity score: " << finalsim << endl;
    return startIndex;
}

void identifyMutations(string input_strand, string target_strand){

    int targetl = target_strand.length();
    int inputl = input_strand.length();

    if(targetl == inputl){
        bool found_mutation = false;

        int index = bestStrandMatch(input_strand, target_strand);
        cout << "Best alignment index: " << index << endl;

        for (int i = 0; i < targetl; i++){
            if(input_strand.at(i) != target_strand.at(i)){
                cout << "Substitution at position " << i + 1 << ": " << input_strand.at(i) << " -> " << target_strand.at(i) << endl;
                found_mutation = true;
            }
        }
        if(found_mutation == false)
        {
            cout << "No mutations found." << endl;
        }
    }

    else if (input_strand.length() < target_strand.length()){

        int index = bestStrandMatch(target_strand, input_strand);
        cout << "Best alignment index: " << index << endl;

        string first = target_strand.substr(0, index);
        string middle = target_strand.substr(index, inputl);
        string last = target_strand.substr(index + inputl);

        if (index > 0){

            for (int i = 0; i < (int)first.length(); i++){
                cout << "Insertion at position " << i + 1 << ": " << first.at(i) << " is inserted in target strand" << endl;
            }

            for (int i = 0; i < (int)middle.length(); i++){
                    if(input_strand.at(i) != middle.at(i)){
                        cout << "Substitution at position " << i + index + 1 << ": " << input_strand.at(i) << " -> " << middle.at(i) << endl;
                    }
            }

            for (int i = 0; i < (int)last.length(); i++){
                cout << "Insertion at position " << i + (int)input_strand.length() + index + 1 << ": " << last.at(i) << " is inserted in target strand" << endl;
            }
        }
        else if(index == 0){

           for (int i = 0; i < (int)middle.length(); i++){
                    if(input_strand.at(i) != middle.at(i)){
                        cout << "Substitution at position " << i + index + 1 << ": " << input_strand.at(i) << " -> " << middle.at(i) << endl;
                    }
            }

            for (int i = 0; i < (int)last.length(); i++){
                cout << "Insertion at position " << i + (int)input_strand.length() + index + 1 << ": " << last.at(i) << " is inserted in target strand" << endl;
            }
        }
    }

    else if(input_strand.length() > target_strand.length()){

        int index = bestStrandMatch(input_strand, target_strand);
        cout << "Best alignment index: " << index << endl;

        string first = input_strand.substr(0, index);
        string middle = input_strand.substr(index, targetl);
        string last = input_strand.substr(index + targetl);

        if (index > 0){

            for (int i = 0; i < (int)first.length(); i++){
                cout << "Deletion at position " << i + 1 << ": " << first.at(i) << " is deleted in target strand" << endl;
            }

            for (int i = 0; i < (int)middle.length(); i++){
                    if(target_strand.at(i) != middle.at(i)){
                        cout << "Substitution at position " << i + index + 1 << ": " << middle.at(i) << " -> " << target_strand.at(i) << endl;
                    }
            }

            for (int i = 0; i < (int)last.length(); i++){
                cout << "Deletion at position " << i + (int)target_strand.length() + index + 1 << ": " << last.at(i) << " is deleted in target strand" << endl;
            }
        }
        else if(index == 0){

           for (int i = 0; i < (int)middle.length(); i++){
                    if(target_strand.at(i) != middle.at(i)){
                        cout << "Substitution at position " << i + index + 1 << ": " << middle.at(i) << " -> " << target_strand.at(i) << endl;
                    }
            }

            for (int i = 0; i < (int)last.length(); i++){
                cout << "Deletion at position " << i + (int)target_strand.length() + index + 1 << ": " << last.at(i) << " is deleted in target strand" << endl;
            }
        }
    }
}

void transcribeDNAtoRNA(string strand){
    int length = strand.length();

    for (int i = 0; i < length; i++){
        if(strand.at(i) == 'T'){
            strand.at(i) = 'U';
        }
    }
    cout << strand << endl;
}

void reverseComplement(string strand){
    int length = strand.length();
    string newstrand = strand;

    for (int i = 0; i < length; i++){
        if (strand.at(i) == 'A'){
            strand.at(i) = 'T';
        }
        else if (strand.at(i) == 'T'){
            strand.at(i) = 'A';
        }
        else if (strand.at(i) == 'C'){
            strand.at(i) = 'G';
        }
        else if (strand.at(i) == 'G'){
            strand.at(i) = 'C';
        }
    }

    for (int i = 0; i < length; i++){

        newstrand.at(i) = strand.at(length - 1 - i);
    }

    cout << newstrand << endl;
}

void getCodingFrames(string strand){

    string newstrand = strand; 
    bool exists = false;

    for (int i = 0; i < (int)strand.length() - 2; i++){
        if(strand.substr(i, 3) == "ATG"){
            newstrand = strand.substr(i, 3);
            i += 3;

            for(int j = i; j < (int)strand.length() - 2; j += 3){

                if(strand.substr(j, 3) == "TAA" || strand.substr(j, 3) == "TAG" || strand.substr(j, 3) == "TGA"){
                    cout << newstrand + strand.substr(i, j + 3 - i) << endl;
                    exists = true;
                    i = j + 2;
                    break;
                }
            }
        }
    }
    if (exists == false){
        cout << "No reading frames found." << endl;
    }
}

bool askPrompt(int choice){

    if(choice == 1 || choice == 2 || choice == 3 || choice == 4 || choice == 5 || choice == 6 || choice == 7){
        return true;
    }
    else{
        cout << "Invalid input. Please select a valid option." << endl;

        cout << "--- DNA Analysis Menu ---" << endl;

        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        return false;
    }
}

void choicenum(int choice){
    if (choice == 1){
        string firstDNA;
        string secondDNA;

        cout << "Enter the first DNA sequence:" << endl;
        cin >> firstDNA;

        while(!isValidStrand(firstDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the first DNA sequence:" << endl;
            cin >> firstDNA;
        }

        cout << "Enter the second DNA sequence:" << endl;
        cin >> secondDNA;

        while(!isValidStrand(secondDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the second DNA sequence:" << endl;
            cin >> secondDNA;
        }

        if((int)firstDNA.length() != (int)secondDNA.length()){
            cout << "Error: Input strands must be of the same length." << endl;
            cout << "--- DNA Analysis Menu ---" << endl;
            cout << "1. Calculate the similarity between two sequences of the same length" << endl;
            cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
            cout << "3. Identify mutations" << endl;
            cout << "4. Transcribe DNA to RNA" << endl;
            cout << "5. Find the reverse complement of a DNA sequence" << endl;
            cout << "6. Extract coding frames" << endl;
            cout << "7. Exit" << endl;
            cout << "Please enter your choice (1 - 7):" << endl;
            cin >> choice;

            while (askPrompt(choice) == false){
                cin >> choice;
                askPrompt(choice);
            }
            if(askPrompt(choice) == true){
                choicenum(choice);
            }

        }
        else{
            cout << "Similarity score: " << strandSimilarity(firstDNA, secondDNA) << endl;
            cout << "--- DNA Analysis Menu ---" << endl;
            cout << "1. Calculate the similarity between two sequences of the same length" << endl;
            cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
            cout << "3. Identify mutations" << endl;
            cout << "4. Transcribe DNA to RNA" << endl;
            cout << "5. Find the reverse complement of a DNA sequence" << endl;
            cout << "6. Extract coding frames" << endl;
            cout << "7. Exit" << endl;
            cout << "Please enter your choice (1 - 7):" << endl;
            cin >> choice;

            while (askPrompt(choice) == false){
                cin >> choice;
                askPrompt(choice);
            }
            if(askPrompt(choice) == true){
                choicenum(choice);
            }
        }
    }

    else if(choice == 2){
        string firstDNA;
        string secondDNA;

        cout << "Enter the first DNA sequence:" << endl;
        cin >> firstDNA;

        while(!isValidStrand(firstDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the first DNA sequence:" << endl;
            cin >> firstDNA;
        }

        cout << "Enter the second DNA sequence:" << endl;
        cin >> secondDNA;

        while(!isValidStrand(secondDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the second DNA sequence:" << endl;
            cin >> secondDNA;
        }

        if ((int)firstDNA.length() < (int)secondDNA.length()){
            bestStrandMatch(secondDNA, firstDNA);
        }
        else{
            bestStrandMatch(firstDNA, secondDNA);
        }

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
        
    }

    else if(choice == 3){
        string firstDNA;
        string secondDNA;

        cout << "Enter the first DNA sequence:" << endl;
        cin >> firstDNA;

        while(!isValidStrand(firstDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the first DNA sequence:" << endl;
            cin >> firstDNA;
        }

        cout << "Enter the second DNA sequence:" << endl;
        cin >> secondDNA;

        while(!isValidStrand(secondDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the second DNA sequence:" << endl;
            cin >> secondDNA;
        }

        identifyMutations(firstDNA, secondDNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
    }
    else if(choice == 4){
        string DNA;

        cout << "Enter the DNA sequence to be transcribed:" << endl;
        cin >> DNA;

        while(!isValidStrand(DNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the DNA sequence to be transcribed:" << endl;
            cin >> DNA;
        }

        cout << "The transcribed DNA is: ";
        transcribeDNAtoRNA(DNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
    }
    else if(choice == 5){
        string DNA;

        cout << "Enter the DNA sequence: " << endl;
        cin >> DNA;

        while(!isValidStrand(DNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the DNA sequence:" << endl;
            cin >> DNA;
        }

        cout << "The reverse complement is: ";
        reverseComplement(DNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
    }
    else if(choice == 6){
        string DNA;

        cout << "Enter the DNA sequence: " << endl;
        cin >> DNA;

        while(!isValidStrand(DNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the DNA sequence:" << endl;
            cin >> DNA;
        }

        cout << "The extracted reading frames are: ";
        getCodingFrames(DNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
    }
    else if(choice == 7){
        cout << "Exiting program." << endl;
    }
}

int main(){

    int choice  = 0;
    
    cout << "--- DNA Analysis Menu ---" << endl;
    cout << "1. Calculate the similarity between two sequences of the same length" << endl;
    cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
    cout << "3. Identify mutations" << endl;
    cout << "4. Transcribe DNA to RNA" << endl;
    cout << "5. Find the reverse complement of a DNA sequence" << endl;
    cout << "6. Extract coding frames" << endl;
    cout << "7. Exit" << endl;
    cout << "Please enter your choice (1 - 7):" << endl;
    cin >> choice;

    while (askPrompt(choice) == false){
        cin >> choice;
        askPrompt(choice);
    }
    

    if (choice == 1){
        string firstDNA;
        string secondDNA;

        cout << "Enter the first DNA sequence:" << endl;
        cin >> firstDNA;

        while(!isValidStrand(firstDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the first DNA sequence:" << endl;
            cin >> firstDNA;
        }

        cout << "Enter the second DNA sequence:" << endl;
        cin >> secondDNA;

        while(!isValidStrand(secondDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the second DNA sequence:" << endl;
            cin >> secondDNA;
        }

        if((int)firstDNA.length() != (int)secondDNA.length()){
            cout << "Error: Input strands must be of the same length." << endl;
            cout << "--- DNA Analysis Menu ---" << endl;
            cout << "1. Calculate the similarity between two sequences of the same length" << endl;
            cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
            cout << "3. Identify mutations" << endl;
            cout << "4. Transcribe DNA to RNA" << endl;
            cout << "5. Find the reverse complement of a DNA sequence" << endl;
            cout << "6. Extract coding frames" << endl;
            cout << "7. Exit" << endl;
            cout << "Please enter your choice (1 - 7):" << endl;
            cin >> choice;

            while (askPrompt(choice) == false){
                cin >> choice;
                askPrompt(choice);
            }
            if(askPrompt(choice) == true){
                choicenum(choice);
            }

        }
        else{
            cout << "Similarity score: " << strandSimilarity(firstDNA, secondDNA) << endl;
            cout << "--- DNA Analysis Menu ---" << endl;
            cout << "1. Calculate the similarity between two sequences of the same length" << endl;
            cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
            cout << "3. Identify mutations" << endl;
            cout << "4. Transcribe DNA to RNA" << endl;
            cout << "5. Find the reverse complement of a DNA sequence" << endl;
            cout << "6. Extract coding frames" << endl;
            cout << "7. Exit" << endl;
            cout << "Please enter your choice (1 - 7):" << endl;
            cin >> choice;

            while (askPrompt(choice) == false){
                cin >> choice;
                askPrompt(choice);
            }
            if(askPrompt(choice) == true){
                choicenum(choice);
            }
        }
    }

    else if(choice == 2){
        string firstDNA;
        string secondDNA;

        cout << "Enter the first DNA sequence:" << endl;
        cin >> firstDNA;

        while(!isValidStrand(firstDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the first DNA sequence:" << endl;
            cin >> firstDNA;
        }

        cout << "Enter the second DNA sequence:" << endl;
        cin >> secondDNA;

        while(!isValidStrand(secondDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the second DNA sequence:" << endl;
            cin >> secondDNA;
        }

        bestStrandMatch(firstDNA, secondDNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
        
    }

    else if(choice == 3){
        string firstDNA;
        string secondDNA;

        cout << "Enter the first DNA sequence:" << endl;
        cin >> firstDNA;

        while(!isValidStrand(firstDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the first DNA sequence:" << endl;
            cin >> firstDNA;
        }

        cout << "Enter the second DNA sequence:" << endl;
        cin >> secondDNA;

        while(!isValidStrand(secondDNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the second DNA sequence:" << endl;
            cin >> secondDNA;
        }

        identifyMutations(firstDNA, secondDNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
    }
    else if(choice == 4){
        string DNA;

        cout << "Enter the DNA sequence to be transcribed:" << endl;
        cin >> DNA;

        while(!isValidStrand(DNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the DNA sequence to be transcribed:" << endl;
            cin >> DNA;
        }

        cout << "The transcribed DNA is: ";
        transcribeDNAtoRNA(DNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
    }
    else if(choice == 5){
        string DNA;

        cout << "Enter the DNA sequence: " << endl;
        cin >> DNA;

        while(!isValidStrand(DNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the DNA sequence:" << endl;
            cin >> DNA;
        }

        cout << "The reverse complement is: ";
        reverseComplement(DNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
    }

    else if(choice == 6){
        string DNA;

        cout << "Enter the DNA sequence: " << endl;
        cin >> DNA;

        while(!isValidStrand(DNA)){
            cout << "Invalid input. Please enter a valid sequence." << endl;
            cout << "Enter the DNA sequence:" << endl;
            cin >> DNA;
        }

        cout << "The extracted reading frames are: ";
        getCodingFrames(DNA);

        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;
        cin >> choice;

        while (askPrompt(choice) == false){
            cin >> choice;
            askPrompt(choice);
        }
        if(askPrompt(choice) == true){
            choicenum(choice);
        }
    }
    else if(choice == 7){
        cout << "Exiting program." << endl;
    }

    return 0;
}