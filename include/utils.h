#pragma once

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Utils (alot of these could be more generalized...)

// Used to Sort those Pesky Reverse Strands
bool compare_first_element(const std::vector<int>& a, const std::vector<int>& b);

// Sort Using Vector Lengths
bool compare_length(const std::vector<int>& a, const std::vector<int>& b);

// Reverse and Negative if reverse strand (I'm actually pretty proud of this solution)
std::vector<int> reverse_and_negate(const std::vector<int> &vec);

// Reverse Transcripts
void reverse_transcripts(std::vector<std::vector<int>> &transcripts);

// Get Min Position of Identified Cluster
int get_pos_min(const int &index, std::vector<std::vector<int>> &core, std::vector<int> *vec);

// Get Max Position of Identified Cluster
int get_pos_max(const int &index, std::vector<std::vector<int>> &core, std::vector<int> *vec);

// Swap Int Variables
void variable_swap(int &a, int &b);

// Checks to see if regions of transcript overlap a gene, also used to overlap transrcipts
bool check_bounds(const int &a_start, const int &a_stop, const int &b_start, const int &b_stop);

// Check Overlap of read position with exon
bool check_point_overlap(const int &p, const int &e1, const int &e2);

// Check if vector is contained within another vector
bool check_containment(const std::vector<int> &b, const std::vector<int> &a);

// For Debugging
void print_transcripts(const std::vector<std::vector<int>> &transcripts);