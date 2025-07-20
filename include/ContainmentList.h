#pragma once

/*
Inspitred by this paper, https://academic.oup.com/bioinformatics/article/23/11/1386/199545
	Although... this technically is not a nested containment list.
*/

class ContainmentList {

public:

	int indices = 0;                     // Number of Starts and Stops
	int regions = 0;                     // Number of bounds (regions of start/stops)
	std::vector<int> vals;               // Bounds

	ContainmentList* next = nullptr;     // Next 
	ContainmentList* prev = nullptr;

	size_t sublist_count = 0;
	std::vector<ContainmentList*> sublist;

	/////////////////////////////////////////////////////////////
	/* Constructors */

	ContainmentList(const std::vector<int> &intervals) {
		for (int i = 0; i < intervals.size(); i++) {
			this -> vals.push_back(intervals[i]);
			++indices;
		}
		regions = indices / 2;
	}

	// Empty
	ContainmentList() {};

	/////////////////////////////////////////////////////////////
	/* List Operations */

	void clean() {
		if (sublist.empty()) { return; }
		for (auto &s : sublist) {
			s -> clean(); delete s;
		}
		sublist.clear();
		sublist_count = 0;
	}

	void set_next(ContainmentList *next) { this -> next = next; }
	int get_back() {
		if (indices == 0) { return -1; }
		return vals[indices - 1];
	}

	void add_interval(const std::vector<int> &intervals) {
		sublist.emplace_back(new ContainmentList(intervals));
		++sublist_count;
	}
	

	void print_intervals() {
		for (int i = 0; i < indices; i++) { std::cerr << vals[i] << ","; }
		std::cerr << "\n";
	}

	/////////////////////////////////////////////////////////////
	/* Transcript Functions */

	std::vector<std::pair<int, int>> make_pairs(const std::vector<int> &a, const std::vector<int> &b);
	std::vector<int> merge_intervals(std::vector<std::pair<int, int>> &pairs);
	void collapse_intervals();
};
