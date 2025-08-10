#include <iostream>
#include <vector>
#include <algorithm>

#include "ContainmentList.h"
#include "global_args.h"
#include "utils.h"

/////////////////////////////////////////////////////////////
/* Transcript Functions */

// Convert two vectors into vector of pairs
std::vector<std::pair<int, int>> ContainmentList::make_pairs(const std::vector<int> &a, const std::vector<int> &b) {
	std::vector<std::pair<int, int>> pairs;
	for (int i = 0; i < a.size() / 2; i++) { pairs.emplace_back(a[(2*i)], a[(2*i)+1]); }
	for (int i = 0; i < b.size() / 2; i++) { pairs.emplace_back(b[(2*i)], b[(2*i)+1]); }
	return pairs;
}


// Merge vector of pairs into nonoverlapping intervals
std::vector<int> ContainmentList::merge_intervals(std::vector<std::pair<int, int>> &pairs) {

	// Sort vector by starts
  std::sort(pairs.begin(), pairs.end(),
            [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                return a.first < b.first;
            });

  // Flatten and merge
  int end;
	std::vector<int> tmp;
  for (int i = 0; i < pairs.size(); i++) {
  	tmp.push_back(pairs[i].first);
  	end = pairs[i].second;
  	
  	for (int j = i + 1; j < pairs.size(); j++) {
  		if (end >= pairs[j].first || std::abs(end - pairs[j].first) <= ImpaqtArguments::Args.epsilon) {
    		end = std::max(end, pairs[j].second);
    		i = j; 
  		} else { break; }
  	}
  	tmp.push_back(end);
  }

	return tmp;
}

// Create Transcript offspring and collapse if possible
void ContainmentList::collapse_intervals() {

	/*
	  The trick here is treat each sub-ContainmentList as it's own potential path.
		We construct the transcript resulting from the parent and each child 
		and then try to merge the children. If children cannot be merged they
		are separated from this family of transcripts and become their own 
		parent in the next evolution of our ContainmentList.
	*/
	
	if (sublist_count == 0) { return; }

	std::vector<std::vector<int>> tmp;
  std::vector<std::pair<int, int>> pairs;

  // Convert to vector of pairs and create transcript
  for (const auto &s : sublist) {
  	pairs = make_pairs(this -> vals, s -> vals);
  	tmp.emplace_back(merge_intervals(pairs));
	}

	// Try to Merge Children
	bool unique;
	while (!unique) {
		
		unique = true;
		for (int i = 0; i < tmp.size(); i++) {
			
			if (tmp[i].empty()) { continue; }

			for (int j = i + 1; j < tmp.size(); j++) {
			
				if (check_containment(tmp[j], tmp[i])) {
					pairs = make_pairs(tmp[i], tmp[j]);
					tmp[i] = merge_intervals(pairs);
					tmp[j] = {};
					unique = false;
				}
			}
		}
	}

	// Create new offspring from merged transcripts.
	this -> clean();
	for (const auto &t : tmp) {
		if (t.empty()) { continue; }
		this -> add_interval(t);
	}
}
