#include <string> // for std::string
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <cmath>
#include <math.h>
#include <vector>
#include <cstdlib> // exit, EXIT_FAILURE
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <set>

#include "fastcluster.lib.cpp"

typedef std::vector<t_float> Values;
typedef std::vector<Values> Matrix;

typedef std::vector<std::string> Labels;
typedef std::vector<size_t> Indices;
typedef std::vector<unsigned> Counts;
typedef std::vector<int> Ints;
typedef std::vector<Ints> Members;
typedef std::vector<Ints> Merge;

static t_float const NANVALUE = std::numeric_limits<t_float>::quiet_NaN();

struct HclustResult {
	Labels labels;
	Merge merge;
	Values height;
	Ints order;
};
typedef std::vector<HclustResult> HclustResults;

std::ostream & operator << (std::ostream & out, Values const & vals ){
	for(Values::const_iterator it(vals.begin()); it!=vals.end(); ++it) out << *it << " ";
	return out;
}

// colinds allows computing pearson over arbitrary combinations of columns (useful for bootstrapping)
void pearson_distances(Matrix const & matrix, Indices const & colinds, t_float * const dist, bool pairwise_complete_obs=true){
	// pairwise pearson correlations
	size_t const nrow(matrix.size()), ncol(colinds.size());
	size_t r1(0), r2(0), i(0);
	std::ptrdiff_t p(0);
	t_float EX(0.0), EY(0.0), EXX(0.0), EYY(0.0), EXY(0.0), x(0.0), y(0.0), d(0.0), denom(0.0);
	unsigned npairs(0);
	for(r1=0; r1<(nrow-1); ++r1){
		for(r2=r1+1; r2<nrow; ++r2){
			EX=0, EY=0, EXX=0, EYY=0, EXY=0;
			npairs=0;
			for(i=0; i<ncol; ++i){
				x = matrix[r1][colinds[i]];
				y = matrix[r2][colinds[i]];
				if(pairwise_complete_obs && (fc_isnan(x) || fc_isnan(y))) continue;
				++npairs;
				EX += x;
				EY += y;
				EXX += x*x;
				EYY += y*y;
				EXY += x*y;
			}

			d = 2.0;
			if(npairs>0){
				// nan values can appear during bootstrap resampling, due to repeated values and a zero denominator
				//d = 1.0 - (EXY - EX*EY/npairs) / sqrt( (EXX - EX*EX/npairs)*(EYY - EY*EY/npairs) );
				denom = (EXX - EX*EX/npairs)*(EYY - EY*EY/npairs);
				if(denom>0) d = 1.0 - (EXY - EX*EY/npairs) / sqrt(denom);
			}

			//std::cerr << d << std::endl;
            if(fc_isnan(d)){
                t_float num(EXY - EX*EY/npairs);
                std::cerr << "nan distance value for r1 " << r1 << " r2 " << r2 << " EX " << EX << " EY " << EY << " EXX " << EXX << " EYY " << EYY << " EXY " << EXY << " npairs " << npairs << " num " << num << " denom " << denom << " sqrtdenom " << sqrt(denom) << std::endl;
                std::cerr << "colinds:";
                for(i=0; i<ncol; ++i) std::cerr << " " << colinds[i];
                std::cerr << std::endl;
            }
			dist[p++] = d;
		}
	}
}

// begin Spearman code
struct OrdVal{
	t_float val;
	size_t orig_ind;
};
bool sortbyval(OrdVal const & ov1, OrdVal const & ov2){
    if(fc_isnan(ov1.val) && fc_isnan(ov2.val)) return false;
    if(fc_isnan(ov1.val) && !fc_isnan(ov2.val)) return false;
    if(!fc_isnan(ov1.val) && fc_isnan(ov2.val)) return true;
    return ov1.val < ov2.val;
}
typedef std::vector<OrdVal> OrdVals;

typedef std::vector<Indices> IndMatrix;

struct Tie{
	size_t index;
	unsigned rank;
};

typedef std::set<unsigned> Ties;

// this function is much, much faster than R's cor(...,method='spearman') because:
// 1) it's c++, 2) it precomputes all ranks prior to N^2 loop.
// handles missing data ('NAs') in a pairwise.complete.obs fashion
void spearman_distances(Matrix const & matrix, Indices const & colinds, t_float * const dist){

	// precompute ranks for all rows
	size_t const nrow(matrix.size()), ncol(colinds.size());
	Matrix ranks(nrow);
	for(size_t row(0); row<nrow; ++row){
        //std::cout << row << std::endl;

		OrdVals ordvals;
		unsigned colcounter(0);
		for(Indices::const_iterator col(colinds.begin()); col!=colinds.end(); ++col){
			t_float val(matrix[row][*col]);
            //std::cout << '\t' << *col << " " << val << std::endl;
			OrdVal ov;
			ov.val = val;
			// don't store *col here, store the original order in which the sampled columns were added
			ov.orig_ind = colcounter;
			ordvals.push_back(ov);
			++colcounter;
		}
        //std::cout << std::endl;
		std::sort(ordvals.begin(), ordvals.end(), sortbyval);

		Indices rankorder;
		for(size_t i(0); i<ordvals.size(); ++i) rankorder.push_back(ordvals[i].orig_ind);

        //for(size_t i(0); i<ordvals.size(); ++i) std::cout << '\t' << i << " " << ordvals[i].val << " " << ordvals[i].orig_ind << std::endl;

		// compute ranks in ascending order (nan's last)
		Values rowranks(ncol,NANVALUE);
		Ties ties;
		t_float last(0);
        unsigned naiverank(0);
		for(size_t i(0); i<ncol; ++i){
            t_float val(ordvals[i].val);
			if(fc_isnan(val)){
                // sort function ensures nan values are at the end, so these can just pass through
				rowranks[i] = val;
				continue;
			}
			// ranks are 1-indexed (though this shouldn't matter?)
			else rowranks[i] = (t_float)(naiverank+1);

			// deal with ties
			if(val == last){
				ties.insert(naiverank+1);
				if(i>0) ties.insert(naiverank);
			} else {
				t_float rank(0);
				for(Ties::const_iterator t(ties.begin()); t!=ties.end(); ++t) rank += *t;
				for(size_t ip(0); ip<ties.size(); ++ip) rowranks[i-ip-1] = rank/ties.size();
				ties.clear();
			}

			last = val;
			++naiverank;
		}

		// process final ties
		if(!ties.empty()){
			t_float rank(0);
			for(Ties::const_iterator t(ties.begin()); t!=ties.end(); ++t) rank += *t;
			for(size_t ip(0); ip<ties.size(); ++ip) rowranks[naiverank-ip-1] = rank / ties.size();
		}

		// now 'unsort' the ranks so that any pair of two ranks will re-correspond to each other
		ranks[row].resize(rankorder.size());
		for(size_t i(0); i<rankorder.size(); ++i) ranks[row][rankorder[i]] = rowranks[i];

        //std::cout << "Ranks:";
        //for(size_t i(0); i<ranks[row].size(); ++i) std::cout << '\t' << i << " " << ranks[row][i] << std::endl;
	}

	// pairwise pearson correlations of the pre-computed spearman ranks
	std::ptrdiff_t p(0);
	t_float EX(0.0), EY(0.0), EXX(0.0), EYY(0.0), EXY(0.0), x(0.0), y(0.0), d(0.0), denom(0.0);
	size_t r1(0), r2(0), i(0);
	unsigned npairs(0);
	for(r1=0; r1<(nrow-1); ++r1){
		for(r2=r1+1; r2<nrow; ++r2){

			EX=0, EY=0, EXX=0, EYY=0, EXY=0;
			npairs=0;
			for(i=0; i<ncol; ++i){
				x = ranks[r1][i];
				y = ranks[r2][i];
				if(fc_isnan(x) || fc_isnan(y)) continue;
				++npairs;
				EX += x;
				EY += y;
				EXX += x*x;
				EYY += y*y;
				EXY += x*y;
			}

			d = 2.0;
			if(npairs>0){
				// nan values can appear during bootstrap resampling, due to repeated values and a zero denominator
				//d = 1.0 - (EXY - EX*EY/npairs) / sqrt( (EXX - EX*EX/npairs)*(EYY - EY*EY/npairs) );
				denom = (EXX - EX*EX/npairs)*(EYY - EY*EY/npairs);
				if(denom>0) d = 1.0 - (EXY - EX*EY/npairs) / sqrt(denom);
			}

			//std::cerr << d << std::endl;
			if(fc_isnan(d)){
				t_float num(EXY - EX*EY/npairs);
				std::cerr << "nan distance value for r1 " << r1 << " r2 " << r2 << " EX " << EX << " EY " << EY << " EXX " << EXX << " EYY " << EYY << " EXY " << EXY << " npairs " << npairs << " num " << num << " denom " << denom << " sqrtdenom " << sqrt(denom) << std::endl;
			}

			dist[p++] = d;

		}
	}
}

void read_matrix_file(
                      std::string const & filename,
                      Matrix & rr,
                      Labels & ids
                      ){
	// open input file
	std::ifstream file;
	file.open( filename.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << filename.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cerr << "Reading file " << filename << std::endl;

	// read input file
	std::string line;
	bool firstline(true);
	while ( getline( file, line ) ) {
		// skip the first line, which should be a header
		if(firstline){ firstline=false; continue; }
		std::istringstream linestream(line);
		Values ratios;
		std::string token;
		bool id(true);
		while(getline(linestream, token, ' ')){
			//			std::cerr << token << ":";
			// first token in line should be an id
			if(id){
				ids.push_back(token);
				id=false;
				continue;
			}
			// the rest should be numbers, but sometimes "NA"
			t_float value;
			if(token=="NA") value=NANVALUE;
			else value=std::atof(token.c_str());
			ratios.push_back(value);
		}

		rr.push_back(ratios);
	}
}

void read_merge_file(std::string const & filename, Merge & merge){
	// open input file
	std::ifstream file;
	file.open( filename.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << filename.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cerr << "Reading file " << filename << std::endl;

	// read input file
	std::string line;
	while (getline(file, line)) {
		std::istringstream linestream(line);
		Ints row;
		int ind;
		while(linestream >> ind) row.push_back(ind);
		merge.push_back(row);
	}
}

// construct a list of hierarchical tree members based on an hclust tree structure (a.k.a. 'merge')
// this is a re-implementation of the hc2split function from pvclust
void get_members(Merge const & merge, Members & members){
	size_t const nrow(merge.size());
	members.clear();
	members.resize(nrow);
	for(size_t i(0); i<nrow; ++i){
		//std::cout << i << " ";
		int ai(merge[i][0]);
		//std::cout << ai << " ";
		if(ai<0) {members[i].clear(); members[i].push_back(-1*ai);}
		// note: indices in R and fastcluster's merge table are 1-indexed
		else members[i] = members[ai-1];

		ai = merge[i][1];
		if(ai<0) members[i].push_back(-1*ai);
		else members[i].insert(members[i].end(),members[ai-1].begin(),members[ai-1].end());
	}
	for(size_t i(0); i<members.size(); ++i)
		std::sort(members[i].begin(), members[i].end());
}

// for counting matching tree nodes (all children indentical) between two hierarchical tree structures ('merges')
// adds to 'nodecounts' (for calculating hclust bootstrap p-values)
void add_equal_nodes(Merge const & merge1, Merge const & merge2, Counts & nodecounts){
	// merge is a two-column matrix of integers from hclust
	size_t const nrow1(merge1.size());

	Members memb1(merge1.size()), memb2(merge2.size());
	get_members(merge1, memb1);
	get_members(merge2, memb2);

	// assumes unique members and does not re-compare to previously matched ones
	std::vector<bool> matched(nrow1,false);
	for(size_t i1(0); i1<memb1.size(); ++i1){
		size_t const msz1(memb1[i1].size());
		for(size_t i2(0); i2<memb2.size(); ++i2){
			if(matched[i2]) continue;
			size_t const msz2(memb2[i2].size());
			if(msz1 != msz2) continue;
			bool match(true);
			for(size_t j(0); j<msz2; ++j){
				if(memb1[i1][j] != memb2[i2][j]) {match=false; break;}
			}
			if(match){
				//std::cout << "counting node at " << i1 << std::endl;
				nodecounts[i1] += 1;
				matched[i2] = true;
			}
		}
	}
}

// output multiple bootstrap results to one file
void output_bootstrap_results(HclustResults const & bootstraps, std::string prefix="bootstrap."){
	// outputs that can be read into R
	std::ofstream of;

	std::string fname = prefix+"labels";
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		std::copy(bs->labels.begin(), bs->labels.end(), std::ostream_iterator<std::string>(of, " "));
		of << '\n';
	}
	of.close();

	fname = prefix+"merges";
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		for(size_t j(0); j<2; ++j){
			for(size_t i(0); i<bs->merge.size(); ++i) of << bs->merge[i][j] << " ";
		}
		of << '\n';
	}
	of.close();

	fname = prefix+"heights";
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		std::copy(bs->height.begin(), bs->height.end(), std::ostream_iterator<t_float>(of, " "));
		of << '\n';
	}
	of.close();

	fname = prefix+"orders";
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		std::copy(bs->order.begin(), bs->order.end(), std::ostream_iterator<int>(of, " "));
		of << '\n';
	}
	of.close();

}

void output_node_counts(Counts const & nodecounts, std::string prefix=""){
	std::ofstream of;
	std::string fname = prefix+"nodecounts";
	of.open(fname.c_str());
	for(Counts::const_iterator it(nodecounts.begin()); it!=nodecounts.end(); ++it) of << *it << '\n';
	of.close();
}

void output_hclust(HclustResult const & result, std::string prefix=""){
	std::ofstream of;

	std::string fname = prefix+"labels";
	of.open(fname.c_str());
	std::copy(result.labels.begin(), result.labels.end(), std::ostream_iterator<std::string>(of, " "));
	of << '\n';
	of.close();

	fname = prefix+"merge";
	of.open(fname.c_str());
	for(size_t j(0); j<2; ++j){
		for(size_t i(0); i<result.merge.size(); ++i) of << result.merge[i][j] << " ";
	}
	of << '\n';
	of.close();

	fname = prefix+"height";
	of.open(fname.c_str());
	std::copy(result.height.begin(), result.height.end(), std::ostream_iterator<t_float>(of, " "));
	of << '\n';
	of.close();

	fname = prefix+"order";
	of.open(fname.c_str());
	std::copy(result.order.begin(), result.order.end(), std::ostream_iterator<int>(of, " "));
	of << '\n';
	of.close();

}

void run_fastcluster(HclustResult & hclust_result, t_float * dist, int method, bool sqrt_heights=false){
	//std::cout << "Fastcluster..." << std::endl;

	if (method<METHOD_METR_SINGLE || method>METHOD_METR_MEDIAN) {
		std::cerr << "ERROR: unknown clustering method requested" << std::endl;
		exit(EXIT_FAILURE);
	}

	const int N((int)hclust_result.labels.size());
	// 'members': members in each node (for 'clustering in the middle of the tree')
	// here, just setting this whole thing to 1 (all data are for terminal branches)
	auto_array_ptr<t_float> members;
	members.init(N);
	for (int i(0); i<N; ++i) members[i] = 1;
	cluster_result hc(N-1);

	switch (method) {

		case METHOD_METR_SINGLE:
			MST_linkage_core(N, dist, hc);
			break;
		case METHOD_METR_COMPLETE:
			NN_chain_core<METHOD_METR_COMPLETE, t_float>(N, dist, NULL, hc);
			break;
		case METHOD_METR_AVERAGE:
			NN_chain_core<METHOD_METR_AVERAGE, t_float>(N, dist, members, hc);
			break;
		case METHOD_METR_WEIGHTED:
			NN_chain_core<METHOD_METR_WEIGHTED, t_float>(N, dist, NULL, hc);
			break;
		case METHOD_METR_WARD:
			NN_chain_core<METHOD_METR_WARD, t_float>(N, dist, members, hc);
			break;
		case METHOD_METR_CENTROID:
			generic_linkage<METHOD_METR_CENTROID, t_float>(N, dist, members, hc);
			break;
		case METHOD_METR_MEDIAN:
			generic_linkage<METHOD_METR_MEDIAN, t_float>(N, dist, NULL, hc);
			break;
		default:
			std::cerr << "ERROR: unknown clustering method requested" << std::endl;
			exit(EXIT_FAILURE);
	}

	// here init vector sizes and pass array pointers to fastcluster
	hclust_result.height.resize(N-1);
	hclust_result.order.resize(N);

	Ints merge(2*(N-1),0);

	if (method==METHOD_METR_CENTROID ||
	    method==METHOD_METR_MEDIAN)
		generate_R_dendrogram<true>(&merge[0], &hclust_result.height[0], &hclust_result.order[0], hc, N);
	else
		generate_R_dendrogram<false>(&merge[0], &hclust_result.height[0], &hclust_result.order[0], hc, N);

    // methods 4, 5, and 6 (Euclidean methods) expect squared distances, and then the R heights must be square rooted to reflect methods accurate to Ward's original method
    if(sqrt_heights){
        for(Values::iterator it(hclust_result.height.begin()); it!=hclust_result.height.end(); ++it)
            *it = sqrt(*it);
    }

	// convert merge (flat array) into two-column Merge matrix for later convenience/intelligibility
	hclust_result.merge.resize(N-1);
	for(int i(0); i<N-1; ++i){
		hclust_result.merge[i].resize(2);
		// merge array is filled by column
		for(size_t j(0); j<2; ++j) hclust_result.merge[i][j] = merge[i+(N-1)*j];
	}

	members.free();
}

void get_distances(
                   Matrix const & matrix,
                   Indices const & colinds,
                   t_float * const dist,
                   std::string method="pearson",
                   bool square_distances=false
                   ){
	if(method=="pearson") pearson_distances(matrix, colinds, dist);
	else if(method=="spearman") spearman_distances(matrix, colinds, dist);
	else{
		std::cerr << "ERROR: Unsupported distance metric " << method << std::endl;
		exit(EXIT_FAILURE);
	}

    if(square_distances){
        const std::ptrdiff_t NN = static_cast<std::ptrdiff_t>((matrix.size())*(matrix.size()-1)/2);
        for(std::ptrdiff_t p(0); p<NN; ++p) dist[p] *= dist[p];
    }
}

std::string method_str(int method){
    std::string str("unknown");
    if(method==0) return "single";
    if(method==1) return "complete";
    if(method==2) return "average";
    if(method==3) return "weighted";
    if(method==4) return "ward";
    if(method==5) return "centroid";
    if(method==6) return "median";
    return str;
}

////////////////////////////////////////////////////////////////////////////////
void usage_error(){
	std::cerr << "\n"
	<< " -r [file]            : matrix of values to cluster, first row is ignored, first column is item labels\n"
	<< " -b #                 : number of bootstraps\n"
	<< " -m [hclust method]   : clustering method (0. single, 1. complete, 2. average, 3. weighted, 4. ward, 5. centroid, 6. median)\n"
    << " -D2                  : flag: square distances? Default is false, but expected by Ward, centroid and median methods in current version (1.1.13) of fastcluster.\n"
	<< " -d [distance metric] : distance metric (supported: pearson or spearman)\n"
	<< " -v #                 : verbosity (1 or 2)\n"
	<< "example: [executable] -r ratios.tab\n"
	<< "\n";
	exit(EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {

	std::srand((unsigned int)time(NULL));

	std::cout << std::endl;

	std::string ratiosfilename, distmethod;
	unsigned bootstraps(0), verbosity(1);
	int method(2); // average
    bool square_distances(false);
	t_float scale(1.0); // ratio of columns to resample for multiscale bootstrapping

	if (argc < 2) usage_error();

	// parse command line arguments
	for (int i(1); i < argc; ++i) {

		std::string arg(argv[i]);

		if (arg == "-r") {
			if (++i >= argc) usage_error();
			ratiosfilename = argv[i];

		} else if (arg == "-b") {
			if (++i >= argc) usage_error();
			bootstraps = atoi(argv[i]);

		} else if (arg == "-s") {
			if (++i >= argc) usage_error();
			scale = atof(argv[i]);

		} else if (arg == "-m") {
			if (++i >= argc) usage_error();
			method = atoi(argv[i]);
            if(method>3) std::cout << "Squared distances are expected for this method: recommend restarting with the -D2 flag" << std::endl;

		} else if (arg == "-d") {
			if (++i >= argc) usage_error();
			distmethod = argv[i];

		} else if (arg == "-D2") {
			square_distances = true;
            std::cout << "Euclidean method: distances will be squared and heights will be square rooted" << std::endl;

		} else if (arg == "-v") {
			if (++i >= argc) usage_error();
			verbosity = atoi(argv[i]);

		} else if (arg == "-h" || arg == "--help") usage_error();
	}

    if(method<4 && square_distances) std::cout << "-D2 flag is set to true but distance method is not expecting squared distances. Recommend restarting without the -D2 flag" << std::endl;

	// read in matrix
	Matrix rr;
	Labels labels;
	read_matrix_file(ratiosfilename, rr, labels);

	if(verbosity>2){
		size_t ntest(std::min(size_t(5),rr.size()));
		std::cout << "Sample data:" << std::endl;
		for(size_t i(0); i<ntest; ++i){ std::cout << labels[i] << " " << rr[i] << std::endl; }
	}

	const size_t nrow(rr.size()), ncol(rr.front().size());

	// dinstance func expects column indices, here just feed 1:ncol (becomes useful for bootstrap resampling)
	Indices colinds;
	for(size_t i(0); i<rr.front().size(); ++i) colinds.push_back(i);

	// Parameter NN: number of non-redundant, non-self comparisons
	const std::ptrdiff_t NN = static_cast<std::ptrdiff_t>((nrow)*(nrow-1)/2);
	std::cout << "NN is " << NN << std::endl;
	std::cout << "Distance matrix..." << std::endl;
	auto_array_ptr<t_float> dist;
	dist.init(NN);

	// feed whole matrix to distance functions:
	// for Spearman, this way can pre-compute ranks prior to the N^2 loop
	get_distances(rr, colinds, dist, distmethod, square_distances);

	// here: do a non-bootstrapped hclust
	HclustResult result;
	// store labels
	for(size_t i(0); i<labels.size(); ++i) result.labels.push_back(labels[i]);

    std::cout << "Fastcluster with method " << method << " (" << method_str(method) << ")" << std::endl;
	run_fastcluster(result, dist, method, square_distances);

	// output result
	output_hclust(result,"hc.");

	// bootstrap iterations
	HclustResults bootstrap_results;
	Counts nodecounts(result.merge.size(),0);

	std::cout << "Bootstrap iterations..." << std::endl;
	unsigned nsample((unsigned)round(ncol*scale));
	for(unsigned bs(0); bs<bootstraps; ++bs){
		if(verbosity>1) std::cout << "Bootstrap iteration " << bs << std::endl;

		HclustResult bs_result;
		// store ids
		for(size_t i(0); i<labels.size(); ++i) bs_result.labels.push_back(labels[i]);

		// resample columns (e.g. microarray conditions) with replacement
		Indices colinds;
		for(t_index i(0); i<nsample; ++i) colinds.push_back(rand() % ncol);

		// this sort probably not necessary
		//std::sort(colinds.begin(), colinds.end());

		if(verbosity>2){
			std::cout << "Colinds: ";
			for(t_index i(0); i<nsample; ++i) std::cout << colinds[i] << " ";
			std::cout << std::endl;
		}

		// resample from columns and create a new distance matrix
		auto_array_ptr<t_float> dist_resampled;
		dist_resampled.init(NN);

        if(verbosity>1) std::cout << "Distance matrix..." << std::endl;
		get_distances(rr, colinds, dist_resampled, distmethod, square_distances);

        if(verbosity>1) std::cout << "Fastcluster..." << std::endl;
		run_fastcluster(bs_result, dist_resampled, method, square_distances);

        dist_resampled.free();

		//bootstrap_results.push_back(bs_result);

		// pvclust-style equal node counting
		add_equal_nodes(result.merge, bs_result.merge, nodecounts);

		// cache current counts
		// outputs that can be read into R
		std::ofstream of;

		std::string fname = "countscurrent";
		of.open(fname.c_str());
		for(Counts::const_iterator it(nodecounts.begin()); it!=nodecounts.end(); ++it) of << *it << '\n';
		of.close();

		fname = "lastiter";
		of.open(fname.c_str());
		of << bs << '\n';
		of.close();

		if(verbosity<2) std::cout << '.' << std::flush;

	}
	std::cout << std::endl;

	output_bootstrap_results(bootstrap_results);
	output_node_counts(nodecounts);

	// maybe in the future: height-based agglomeration to yield k clusters (see dendextendRcpp?),
	// empirical distributions of pairwise cluster co-memberships for many choices of k
	std::cout << "Finished." << std::endl;
	return 0;
}

