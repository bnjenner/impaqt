//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Utils (alot of these could be more generalized...)

// Used to Sort those Pesky Reverse Strands
bool compare_first_element(const std::vector<int>& a, const std::vector<int>& b) { return a[0] < b[0]; }

// Sort Using Vector Lengths
bool compare_length(const std::vector<int>& a, const std::vector<int>& b) { return a.size() < b.size(); }

// Checks to see if regions of transcript overlap
bool check_bounds(const int &a_start, const int &a_stop, const int &b_start, const int &b_stop) {
	// Precedes
	if (a_stop >= b_start && a_stop <= b_stop) {
		return true;
		// Spans
	} else if (a_start <= b_start && a_stop >= b_stop) {
		return true;
		// Follows
	} else if (a_start >= b_start && a_start <= b_stop) {
		return true;
		// Precedes (but close enough)
	} else if (std::abs(a_stop - b_start) <= ImpaqtArguments::Args.epsilon) {
		return true;
	}
	return false;
}

// Reverse and Negative if reverse strand (I'm actually pretty proud of this solution)
std::vector<int> reverse_and_negate(const std::vector<int> &vec) {
	int x = 0;
	int n = vec.size();
	std::vector<int> tmp_vec(n);
	for (int j = n - 1; j > -1; j--) {
		tmp_vec[x] = vec[j] * -1;
		++x;
	}
	return tmp_vec;
}

// Reverse Transcripts
void reverse_transcripts(std::vector<std::vector<int>> &transcripts) {
	for (int i = 0; i < transcripts.size(); i++) {
		transcripts[i] = reverse_and_negate(transcripts[i]);
	}
	std::sort(transcripts.begin(), transcripts.end(), compare_first_element);
}

// Get Min Position of Identified Cluster
int get_pos_min(const int &index, std::vector<std::vector<int>> &core, std::vector<int> *vec) {
	std::vector<int>::iterator min_result = std::min_element(core.at(index).begin(), core.at(index).end());
	return vec -> at(*min_result);
}

// Get Max Position of Identified Cluster
int get_pos_max(const int &index, std::vector<std::vector<int>> &core, std::vector<int> *vec) {
	std::vector<int>::iterator max_result = std::max_element(core.at(index).begin(), core.at(index).end());
	return vec -> at(*max_result);
}

// Swap Int Variables
void variable_swap(int &a, int &b) {
	if (a > b) {
		b = b + a;
		a = b - a;
		b = b - a;
	}
}