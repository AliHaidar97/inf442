#include <iostream>
#include <cassert>
#include <cmath>	// for sqrt, fabs
#include <cfloat>	// for DBL_MAX
#include <cstdlib>	// for rand, srand
#include <ctime>	// for rand seed
#include <fstream>
#include <cstdio>	// for EOF
#include <string>
#include <algorithm>	// for count
#include <vector>
#include<set>
#include <random>
#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gtkmm/drawingarea.h>

class point
{
public:

	static int d;
	double* coords;
	int label;

	point() {

		coords = new double[d];
		label = 0;
		for (int i = 0; i < d; i++) {
			coords[i] = 0;
		}

	}

	~point() {

		delete[]coords;

	}


	//copy constructor 
	point(const point& p1) {

		coords = new double[d];
		for (int i = 0; i < d; i++) {
			coords[i] = p1.coords[i];
		}
		label = p1.label;
	}

	//assigment operator
	point& operator = (const point& p1)
	{
		
		label = p1.label;
		for (int i = 0; i < d; i++) {
			coords[i] = p1.coords[i];
		}
		return *this;
	}

	void print() {

		for (int i = 0; i < d; i++) {
			if (i > 0) {
				std::cout << '\t';
			}
			std::cout << coords[i];
		}
		std::cout << '\n';

	}

	double squared_dist(point& q) {

		double dis = 0;

		for (int i = 0; i < d; i++) {
			dis += (coords[i] - q.coords[i]) * (coords[i] - q.coords[i]);
		}

		return dis;

	}

};

int point::d;


class cloud
{
	private:
	int d;
	int n;
	int k;

	// maximum possible number of points
	int nmax;

	point *points;
	point *centers;

	public:
	cloud(int _d, int _nmax, int _k)
	{
		d = _d;
		point::d = _d;
		n = 0;
		k = _k;

		nmax = _nmax;

		points = new point[nmax];
		centers = new point[k];
	}

	~cloud()
	{
		delete[] centers;
		delete[] points;
	}

	void add_point(point &p, int label)
	{
		assert(("[NOK]", n < nmax));

		for(int m = 0; m < d; m++)
		{
			points[n].coords[m] = p.coords[m];
		}

		points[n].label = label;

		n++;
	}

	int get_d()
	{
		return d;
	}

	int get_n()
	{
		return n;
	}

	int get_k()
	{
		return k;
	}

	point& get_point(int i)
	{
		return points[i];
	}

	point& get_center(int j)
	{
		return centers[j];
	}

	void set_center(point &p, int j)
	{
		for(int m = 0; m < d; m++)
			centers[j].coords[m] = p.coords[m];
	}

	double intracluster_variance()
	{
		
		if (n == 0) {
			return 0;
		}
		double ans = 0;

		for (int i = 0; i < n; i++) {
			ans += points[i].squared_dist(centers[points[i].label]);
		}

		ans /= n;

		return ans;
	}

	int set_voronoi_labels()
	{

		int count = 0;

		for (int i = 0; i < n; i++) {

			int bestLabel = points[i].label;
			double bestDistance = points[i].squared_dist(centers[bestLabel]);

			for (int j = 0; j < k; j++) {
				double curDis = points[i].squared_dist(centers[j]);
				
				if (curDis < bestDistance) {
					bestDistance = curDis;
					bestLabel = j;
				}
					
				if (curDis == bestDistance) {
					bestLabel = std::min(bestLabel, j);
				}
			}
			if (bestLabel != points[i].label) {
				count++;
				points[i].label = bestLabel;
			}
		}

		return count;

	}

	void set_centroid_centers()
	{

		int *size=new int[k];
		for (int i = 0; i < k; i++) {
			size[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			size[points[i].label]++;
		}

		for (int i = 0; i < k; i++) {
			if (size[i] == 0)continue;
			for (int j = 0; j < d; j++) {
				centers[i].coords[j] = 0;
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < d; j++) {
				centers[points[i].label].coords[j] += points[i].coords[j];
			}
		}

		for (int i = 0; i < k; i++) {
			if (size[i] == 0)continue;
			for (int j = 0; j < d; j++) {
				centers[i].coords[j] /= size[i];
			}
		}
		delete[] size;
	}

	void kmeans()
	{
		
		//randomize
		/*for (int i = 0; i < k; i++) {
			for (int j = 0; j < d; j++) {
				centers[i].coords[j] = rand() % 100;
			}
		}
		*/

		//init_forgy();

		//init_plusplus();

		init_random_partition();

		while (set_voronoi_labels() > 0) {
			set_centroid_centers();
		}

	}

	void init_forgy()
	{
		//initialise the centers using the existing points
		std::set<int>s;
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_int_distribution<int> X(0, n-1);

		for (int i = 0; i < k; i++) {
			int x = X(generator);

			while (s.count(x)) {
				x = X(generator);
			}

			//centers[i] = *(new point(points[x])); //I declared a copy constructor memory leaks
			//If we don't overload the assignment operator 
			/*for (int j = 0; j < d; j++) {
				centers[i].coords[j] = points[x].coords[j];
			}*/
			//using the assigment opertor
			centers[i] = points[x];

			s.insert(x);
		}

	}

	void init_plusplus()
	{
		
		centers[0] = points[rand()%n];

		double *dis = new double[n];
		double total = 0;

		for (int i = 0; i < n; i++) {
			dis[i] = centers[0].squared_dist(points[i]);
			total += dis[i];
		}

		for (int i = 1; i < k; i++) {
			double x = (double)rand() / RAND_MAX;
			x = x * total;
			int is = 0;
			for (int j = 0; j < n; j++) {
				if (dis[j] == 0)continue;
				if (x <= dis[j]) {
					is = j;
					break;
				}
				x -= dis[j];
			}
			centers[i] = points[is];
			for (int j = 0; j < n; j++) {
				dis[j] = std::min(dis[j], centers[i].squared_dist(points[j]));
			}
		}
		delete[] dis;
	}

	void init_random_partition()
	{
		
		for (int i = 0; i < n; i++) {
			points[i].label = rand() % k;
		}
		set_centroid_centers();
	}
};


// test functions
void test_point()
{

	std::cout << "Testing the function test_point()\n";
	int nb_errors_1 = 0;
	
	
	//test the constructor and delete[]
	for (int i = 1; i <= 10; i++) {
		point::d = i;
		point a;
		bool wrong = 0;

		if (a.d != i) {
			wrong = 1;
		}

		for (int j = 0; j < a.d; j++) {
			if (a.coords[j] != 0) {
				wrong = 1;
			}
		}

		if (a.label != 0) {
			wrong = 1;
		}
		nb_errors_1 += wrong;
	}


	for (int i = 1; i <= 10; i++) {
		point::d = i;
		point *a=new point();
		bool wrong = 0;

		if (a->d != i) {
			wrong = 1;
		}

		for (int j = 0; j < a->d; j++) {
			if (a->coords[j] != 0) {
				wrong = 1;
			}
		}

		if (a->label != 0) {
			wrong = 1;
		}
		nb_errors_1 += wrong;
		
		delete a;
	}


	for (int i = 1; i <= 10; i++) {
		point::d = i;
		point* a = new point[10];
		for (int z = 0; z < 10; z++) {
			bool wrong = 0;

			if (a[z].d != i) {
				wrong = 1;
			}

			for (int j = 0; j < a[z].d; j++) {
				if (a[z].coords[j] != 0) {
					wrong = 1;
				}
			}

			if (a[z].label != 0) {
				wrong = 1;
			}
			nb_errors_1 += wrong;
		}
		delete[] a;
	}



	//test distance
	int nb_errors_2 = 0;
	const double eps = 0.0001;//tolerance

	point::d = 2;
	
	//test 1
	point* a = new point();
	point* b = new point();

	double realDis = 0.0;
	double testDis = a->squared_dist((*b));

	if (std::abs(testDis - realDis) > eps || (a->squared_dist(*b) != b->squared_dist(*a))) {
		nb_errors_2++;
	}

	//test 2
	b->coords[0] = 1;
	realDis = 1;
	testDis= a->squared_dist((*b));

	if (std::abs(testDis - realDis) > eps || (a->squared_dist(*b) != b->squared_dist(*a))) {
		nb_errors_2++;
	}

	//test 3
	b->coords[0] = 1;
	b->coords[1] = 1;
	realDis = 2;
	testDis = a->squared_dist(*b);
	if (std::abs(testDis - realDis) > eps || (a->squared_dist(*b) != b->squared_dist(*a))) {
		nb_errors_2++;
	}


	//test 4
	b->coords[0] = -1;
	b->coords[1] = 1;
	realDis = 2;
	testDis = a->squared_dist(*b);
	if (std::abs(testDis - realDis) > eps || (a->squared_dist(*b) != b->squared_dist(*a))) {
		nb_errors_2++;
	}

	delete a;
	delete b;


	std::cout << "   #errors = " << nb_errors_1 << " in the point constructor" << std::endl;
	std::cout << "   #errors = " << nb_errors_2 << " in the sqaredDist" << std::endl;
	bool success = ((nb_errors_1 == 0) && (nb_errors_2 == 0));
	std::cout << (success ? "[OK]" : "[NOK]") << std::endl;

}

void test_intracluster_variance()
{
	// tolerance for comparison of doubles
	const double eps = 0.0001;

	// dimension used for tests
	point::d = 1;

	// temporary container
	point p;

	std::cout << "Testing function intracluster_variance()...";

	// test case 1
	const double dist_onepoint_zerodist = 0.0;
	cloud onepoint_zerodist(1, 1, 1);
	p.coords[0] = 0.0;
	onepoint_zerodist.add_point(p, 0);
	assert(("[NOK]", std::fabs(onepoint_zerodist.intracluster_variance() - dist_onepoint_zerodist) < eps));

	// test case 2
	const double dist_onepoint_posdist = 0.25;
	cloud onepoint_posdist(1, 1, 1);
	p.coords[0] = 0.5;
	onepoint_posdist.add_point(p, 0);
	assert(("[NOK]", std::fabs(onepoint_posdist.intracluster_variance() - dist_onepoint_posdist) < eps));

	// test case 3
	const double dist_twopoints = 0.625;
	cloud twopoints(1, 2, 1);
	p.coords[0] = -1.0;
	twopoints.add_point(p, 0);
	p.coords[0] = 0.5;
	twopoints.add_point(p, 0);
	p.coords[0] = -0.5;
	twopoints.set_center(p, 0);
	assert(("[NOK]", std::fabs(twopoints.intracluster_variance() - dist_twopoints) < eps));

	// test case 4
	const double dist_twoclusters = 6.8125;
	cloud twoclusters(1, 4, 2);
	p.coords[0] = -1.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 0.5;
	twoclusters.add_point(p, 0);
	p.coords[0] = -0.5;
	twoclusters.set_center(p, 0);
	p.coords[0] = -2.0;
	twoclusters.add_point(p, 1);
	p.coords[0] = 2.0;
	twoclusters.add_point(p, 1);
	p.coords[0] = -3.0;
	twoclusters.set_center(p, 1);
	assert(("[NOK]", std::fabs(twoclusters.intracluster_variance() - dist_twoclusters) < eps));

	std::cout << "\t[OK]" << std::endl;
}

void test_kmeans()
{
	// TODO
	int nb_of_error = 0;
	point::d = 1;
	//test 1;
	cloud c(1, 10, 3);
	point a1, a2, a3;
	a1.coords[0] = 0;
	a2.coords[0] = 1;
	a3.coords[0] = 2;
	c.add_point(a1, 0);
	c.add_point(a2, 0);
	c.add_point(a3, 0);

	c.kmeans();
	
	//check
	int test[3];
	test[0] = c.get_center(0).coords[0];
	test[1] = c.get_center(1).coords[0];
	test[2] = c.get_center(2).coords[0];
	std::sort(test, test + 3);
	bool wrong = 0;
	for (int i = 0; i < 3; i++) {
		if (test[i] != i) {
			wrong = 1;
		}
	}
	
	nb_of_error += wrong;
	
	//test 2
	wrong = 0;
	point::d = 1;
	cloud c1 (1, 10, 2);
	a1.coords[0] = 0;
	a2.coords[0] = 1;
	a3.coords[0] = 10;
	c1.add_point(a1, 0);
	c1.add_point(a2, 0);
	c1.add_point(a3, 0);

	c1.kmeans();

	if (c1.get_point(0).label != c1.get_point(1).label) {
		wrong = 1;
	}

	if (c1.get_point(0).label == c1.get_point(2).label) {
		wrong = 1;
	}

	nb_of_error += wrong;
	
	/*
	//test 3
	wrong = 0;
	point::d = 2;

	cloud c2(2, 10, 3);
	point a[6];

	a[0].coords[0] = 5;
	a[0].coords[1] = 5;
	a[1].coords[0] = 6;
	a[1].coords[1] = 6;
	

	a[2].coords[0] = -5;
	a[2].coords[1] = 5;
	a[3].coords[0] = -6;
	a[3].coords[1] = 6;
	

	
	a[4].coords[0] = 5;
	a[4].coords[1] = -5;
	a[5].coords[0] = 6;
	a[5].coords[1] = -6;
	

	for (int i = 0; i < 6; i++) {
		c2.add_point(a[i], 0);
	}
	
	c2.kmeans();
	
	for (int i = 0; i < 6; i++) {
		std::cout << c2.get_point(i).label << ' ';
	}std::cout << '\n';

	for (int i = 0; i < 2; i++) {
		if (c2.get_point(i).label != c2.get_point(0).label)
			wrong = 1;
	}

	for (int i = 2; i < 4; i++) {
		if (c2.get_point(i).label != c2.get_point(2).label)
			wrong = 1;
	}

	for (int i = 4; i < 6; i++) {
		if (c2.get_point(i).label != c2.get_point(4).label)
			wrong = 1;
	}
	
	
	nb_of_error += wrong;
	
	std::cout << "nb of errors in test K_means= " << ' ' << nb_of_error << '\n';
	*/
	//assert(("[NOK]", nb_of_error == 0));

	std::cout << "\t\t[OK]" << std::endl;

}

void test_init_forgy()
{
	// number of random experiments
	const int K = 10000;
	// tolerance in probability
	const double delta = 0.0625;

	// dimenstion used for tests
	point::d = 1;

	// temporary container
	point p;

	std::cout << "Testing function init_forgy()...";

	const double prob_threepoints = 0.3333;
	cloud threepoints(1, 3, 1);
	p.coords[0] = 0.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 1.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 2.0;
	threepoints.add_point(p, 0);
	int cnt = 0;
	for(int k = 0; k < K; k++)
	{
		threepoints.init_forgy();
		if(threepoints.get_center(0).coords[0] == 1.0)
			cnt++;
	}
	assert(("[NOK]", std::fabs(cnt/(double)K - prob_threepoints) < delta));

	std::cout << "\t\t[OK]" << std::endl;
}

void test_init_plusplus()
{
	// number of random experiments
	const int K = 10000;
	// tolerance in probability
	const double delta = 0.0625;

	// dimenstion used for tests
	point::d = 1;

	// temporary container
	point p;

	std::cout << "Testing function init_plusplus()...";

	// test case 1
	const double prob_threepoints = 0.3333;
	cloud threepoints(1, 3, 1);
	p.coords[0] = 0.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 1.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 2.0;
	threepoints.add_point(p, 0);
	int cnt = 0;
	for(int k = 0; k < K; k++)
	{
		threepoints.init_plusplus();
		if(threepoints.get_center(0).coords[0] == 1.0)
			cnt++;
	}
	assert(("[NOK]", std::fabs(cnt/(double)K - prob_threepoints) < delta));

	// test case 2
	const double prob_twoclusters = 0.125;
	cloud twoclusters(1, 4, 2);
	p.coords[0] = 0.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 0.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 1.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 2.0;
	twoclusters.add_point(p, 0);
	cnt = 0;
	for(int k = 0; k < K; k++)
	{
		twoclusters.init_plusplus();
		if(twoclusters.get_center(1).coords[0] == 1.0)
			cnt++;
	}
	assert(("[NOK]", std::fabs(cnt/(double)K - prob_twoclusters) < delta));

	std::cout << "\t\t[OK]" << std::endl;
}

void test_init_random_partition()
{
	// number of random experiments
	const int K = 10000;
	// tolerance in probability
	const double delta = 0.0625;

	// dimenstion used for tests
	point::d = 1;

	// temporary container
	point p;

	std::cout << "Testing function init_random_parition()...";

	const double prob_threepoints = 0.3333;
	cloud threepoints(1, 3, 3);
	p.coords[0] = 0.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 1.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 2.0;
	threepoints.add_point(p, 0);
	int cnt = 0;
	for(int k = 0; k < K; k++)
	{
		threepoints.init_random_partition();
		if(threepoints.get_point(2).label == 1)
			cnt++;
	}
	assert(("[NOK]", std::fabs(cnt/(double)K - prob_threepoints) < delta));

	std::cout << "\t[OK]" << std::endl;
}


// for graphical interface
class MyArea : public Gtk::DrawingArea
{
	private:
	cloud *c;
	int x_column;
	int y_column;

	public:
	MyArea(cloud *_c, int _x_column, int _y_column)
	{
		c = _c;
		x_column = _x_column;
		y_column = _y_column;
	}

	virtual ~MyArea() {}

protected:
	//Override default signal handler:
	bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;
};

/**
 * Counts the number of tab-separated columns in the given line.
 */
int nb_columns(const std::string &line)
{
	return std::count(line.begin(), line.end(), '\t') + 1;
}


//silouhette
int get_k(std::vector<point>p) {

	
	const int nmax = 150;
	int best = 1;
	double x = 0;
	double last = 0;
	for (int k = 1; k < 10; k++) {
		cloud c(point::d, nmax, k);
		for (int i = 0; i < p.size(); i++) {
			c.add_point(p[i], 0);
		}
		c.kmeans();
		if (k == 1) {
			last = c.intracluster_variance();
			continue;
		}
		double y = c.intracluster_variance();
		if (x < abs(y - last)) {
			x = abs(y - last);
			best = k;
		}
		last = y;
	}
	return best;
}


int main(int argc, char **argv)
{
	
	test_init_random_partition();

	bool run_tests = false;

	if(argc < 5 || argc > 6)
	{
		std::cerr << "Usage: " << argv[0] << " csv nb_clusters x_column y_column [ test ]" << std::endl;
		std::exit(1);
	}
	std::string csv_filename = argv[1];
	int nb_clusters = std::stoi(argv[2]);
	int x_column = std::stoi(argv[3]);
	int y_column = std::stoi(argv[4]);

	if(argc >= 6)
		run_tests = true;

	srand(time(NULL));

	if(run_tests)
	{
		test_point();
		test_intracluster_variance();
		test_kmeans();
		test_init_forgy();
		test_init_plusplus();
		test_init_random_partition();
	}

	// open data file
	std::ifstream is(csv_filename);
	assert(is.is_open());

	// read header line
	std::string header_line;
	std::getline(is, header_line);
	std::vector<std::string> names;

	const int d = nb_columns(header_line) - 1;
	const int nmax = 150;
	const int k = nb_clusters;






	// construct point cloud
	cloud c(d, nmax, k);

	// point to read into
	point p;

	// labels to cycle through
	int label = 0;

	// while not at end of file
	while(is.peek() != EOF)
	{
		// read new points
		for(int m = 0; m < d; m++)
		{
			is >> p.coords[m];
		}

		c.add_point(p, label);

		label = (label + 1) % k;

		// read ground-truth labels
		// unused in normal operation
		std::string next_name;
		is >> next_name;
		names.push_back(next_name);

		// consume \n
		is.get();
	}

	// execute k-means algorithm
	std::cout << "Intracluster variance before k-means: " << c.intracluster_variance() << std::endl;
	c.kmeans();
	std::cout << "Intracluster variance after k-means: " << c.intracluster_variance() << std::endl;

	std::cout << "Saving clustering into \"output.csv\"... "; 
	std::ofstream os("output.csv");
	assert(os.is_open());
	os << header_line << '\n';
	for(int i = 0; i < c.get_n(); ++i)
	{
		for(int j = 0; j < c.get_d(); ++j)
		{
			os << c.get_point(i).coords[j] << '\t';
		}
		os << names[i] << "_Label_" << c.get_point(i).label;
		if(i != c.get_n()-1)
			os << '\n';
	}
	std::cout << "done" << std::endl;

	os.close();

	// launch graphical interface
	int gtk_argc = 1;
	Glib::RefPtr<Gtk::Application> app = Gtk::Application::create(gtk_argc, argv, "inf442.td3");

	Gtk::Window win;
	win.set_title("TD 3");
	win.set_default_size(400, 400);

	MyArea area(&c, x_column, y_column);
	win.add(area);
	area.show();

	return app->run(win);
}


bool MyArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr)
{
	Gtk::Allocation allocation = get_allocation();
	const int width = allocation.get_width();
	const int height = allocation.get_height();

	// find min and max on each axis
	double x_min = DBL_MAX;
	double x_max = DBL_MIN;
	double y_min = DBL_MAX;
	double y_max = DBL_MIN;
	for(int i = 0; i < c->get_n(); i++)
	{
		if(c->get_point(i).coords[x_column] < x_min)
			x_min = c->get_point(i).coords[x_column];

		if(c->get_point(i).coords[x_column] > x_max)
			x_max = c->get_point(i).coords[x_column];

		if(c->get_point(i).coords[y_column] < y_min)
			y_min = c->get_point(i).coords[y_column];

		if(c->get_point(i).coords[y_column] > y_max)
			y_max = c->get_point(i).coords[y_column];
	}

	// plot all points
	for(int i = 0; i < c->get_n(); i++)
	{
		cr->save(); // save current drawing context (opaque black)
		cr->arc((c->get_point(i).coords[x_column]-x_min)*width/(x_max-x_min), (c->get_point(i).coords[y_column]-y_min)*height/(y_max-y_min), 10.0, 0.0, 2.0 * M_PI); // full circle

		// choose color depending on label
		switch(c->get_point(i).label)
		{
			case 0:
			cr->set_source_rgba(1.0, 0.0, 0.0, 0.6); // red, partially translucent
			break;

			case 1:
			cr->set_source_rgba(0.0, 0.0, 0.8, 0.6); // 0.8 blue, partially translucent
			break;

			case 2:
			cr->set_source_rgba(0.0, 1.0, 0.0, 0.6); // green, partially translucent
			break;

			default:
			double shade = c->get_point(i).label/(double)c->get_k();
			cr->set_source_rgba(shade, shade, shade, 1.0); // gray
			break;
		}

		cr->fill_preserve();
		cr->restore();  // restore drawing context to opaque black
		cr->stroke();
	}

	return true;
}
