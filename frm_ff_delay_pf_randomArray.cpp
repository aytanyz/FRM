
/********************		Gene expression with delay		**************************/

// g++ -std=c++17 -Wall -g -I ~/fastflow -pthread -DNDEBUG -DTRACE_FASTFLOW -finline-functions -O3 -fopenmp frm_ff_delay.cpp -o ff

#include <iostream>
#include <math.h>
#include <limits.h>
#include <vector>
#include <string>
#include <omp.h>
#include <chrono>
#include <utility> 
#include <time.h> 
#include <mutex>
#include <random>
#include <time.h>
#include <algorithm>
#include <stdlib.h>
#include <thread>
#include <atomic>
#include <pthread.h>

#include <ff/ff.hpp>
#include <ff/pipeline.hpp>
#include <ff/parallel_for.hpp>

#include "utimer.cpp"

using namespace std;
using namespace ff;

#define ENDTIME   2				// end of time
#define M 		  3				// number of reaction (should be divisible by 3)
#define N		  4				// number of species (should be divisible by 4)
#define tDelay	  0.1		// delay time for reactions NCD, CD
#define arraysize 980000





// First Stage. Emitter of the Farm
struct Emitter: ff_node_t< double >
{	
	int num_worker;
	
	Emitter(int n_worker) 
	{
		this->num_worker = n_worker;
	}
	
	
	int svc_init()
	{
		// in order to give each thread the different seed
		auto thread_id = std::this_thread::get_id();
		uint64_t* ptr=(uint64_t*) &thread_id;
		double seed = *ptr + time(NULL);		
		srand(seed);
		return 0;
	}
	
	
	double* svc(double *) 
	{
		//cout<<"Generating random seeds for each simulation"<<endl;
		int ___n=num_worker;
		
		for(int i=0; i<___n; i++)
		{
			double rand_seed = (double)rand();;
			ff_send_out(new double(rand_seed));
		}
			
		return EOS;
	};
};

// Second Stage. Worker of the Farm
struct Simulation: ff_node_t< double , vector<int> > 
{
	int 	num_par_for ;
	int 	num_worker;
	int		num_simulation;
	
	int 	x[M];	 			// population of chemical species
	double 	p[N];		    	// propensities of reactants
	
	double 	rate[N];		    // reaction rates
	int		reactant[N][M];
	int		product[N][M];
	int		withDelay[N];		//if reaaction is CD = 1 , NCD = 2, ND = 3 
	double 	myrandQueue[arraysize];
	
	Simulation(int n_worker, int n_simulation, int n_par_for)
	{
		this->num_worker = n_worker;
		this->num_simulation = n_simulation;
		this->num_par_for = n_par_for;
	}
	
	vector<int>* svc(double *seed)
	{		
		//utimer t0("Worker svc: ");
		double 	&my_seed = *seed;		
		vector<int> 	simulation;
		
		int __nw = num_worker;
		int __ns = num_simulation;
		int __npf = num_par_for;
		
		
		
		for(int k=0; k<M; k++)
			simulation.push_back(0); 
		
		init();	
		int count=0;
		generateRandomArrays(my_seed, __npf);
		
		
		
		/************** Starting  num_simulation/num_worker number of simulations************/
		for(int sim=0; sim<__ns/__nw; sim++)
		{
			auto aSim = chrono::high_resolution_clock::now();
			auto eRand =0;
			queue  	<pair<double, int>> q; 				// this pair saves excactly at what time  reaction will finish, and the id of this reaction 
			double 	t =0.0;
		
			x[0] = 1000;
			x[1] = 0;	
			x[2] = 0;
			//auto aIter = chrono::high_resolution_clock::now();
			//cout<<"Starting simulation: "<<sim<<endl;		
			while(t<ENDTIME)
			{	
				//cout<<"sim "<<sim<<" t="<<t<<endl;
				double 	p_sum = propensity_sum();				
				double	min_val = INT_MAX;
				int		min_val_id = -1;
				
				auto aRand = chrono::high_resolution_clock::now();
				/************** Generating N random numbers and finfing smallest one ************/
				for(int i=0; i<N; i++)
				{
					if(count==arraysize)
					{
						//cout<<"count ="<<count<<endl;
						count=0;
						generateRandomArrays(my_seed, __npf);
					}
					
					double myrand = myrandQueue[count];
					count++;
					//cout<<"my rand ="<<myrandQueue[count]<<endl;
					//cout<<"my rand: "<<myrand<<endl; 
					//cout.precision(20);
					//cout<<"sim: "<<sim<<"& id :"<<id<<" and rand: "<<myrand<<endl;
					double randm = -log(myrand) / p_sum;
					//cout<<"random: "<<randm<< endl;
					
					//cout.precision(20);
					//cout<<"Generated pair: "<<randm<<" , "<<i<<endl;
					
					if(randm<min_val)
					{
						min_val = randm;
						min_val_id = i;
					}
				}
				auto myT = chrono::high_resolution_clock::now() - aRand;
				eRand +=  chrono::duration_cast<chrono::nanoseconds>(myT).count();
				 
				
				double	timeTau 	= min_val;
				int 	idOfReac	= min_val_id;
				bool 	chechToFire = true;
				
					
				/************** Checking reactions on the waiting queue ************/
				bool reactionFiredFromQueue = false;
				if(!q.empty())
				{
					//cout<<"Queue is not empty!"<<endl;
					
					std::pair<double, int> previousReaction = q.front();
					
					if(previousReaction.first <= t+timeTau)
					{
						//cout<<"Reaction from queue can be fired now."<<endl;
						
						//update propensities of previousReaction
						switch(withDelay[previousReaction.second])
						{
							case 1: //CD	
										{
											// update the abundance of products
											for(int i=0; i<M; i++)
												x[i] += product[previousReaction.second][i]; 
											break;
										}
							case  2: //NCD	
										{
											// update the abundance of reactants and products
											for(int i=0; i<M; i++)
											{
												x[i] += reactant[previousReaction.second][i] + product[previousReaction.second][i];
											}
											break;
										}
						}
						
						t = previousReaction.first;
						q.pop();
						reactionFiredFromQueue = true;
					}
				}
						
				/************** Firing the chosen reaction ************/
				//cout<<"Now let see if we can fire our new reaction ......"<<endl;
				
				// First we check if there are enough reactants to fire this reaction
				chechToFire = Check_reaction(idOfReac);
				
				if(reactionFiredFromQueue == false && chechToFire) // propensities are not changed, so we can fire our new reaction
				{
					switch(withDelay[idOfReac])
					{
						case 1: //CD 
									{
										// update the abundance of reactants
										for(int i=0; i<M; i++)
											x[i] += reactant[idOfReac][i];	
										//add to the queue
										q.push(make_pair(t+tDelay, idOfReac));										
										//cout<<"CD reaction goes to queue becaus of delay."<<endl;
										break;
									}
						case 2: //NCD	
									{	
										//add to the queue
										q.push(make_pair(t+tDelay, idOfReac));										
										//cout<<"NCD reaction goes to queue becaus of delay."<<endl;
										break;
									}
						case 3: //ND	
									{
										// update the amount of species
										for(int i=0; i<M; i++)
											x[i] += reactant[idOfReac][i] + product[idOfReac][i];										
										//cout<<"ND reaction fires."<<endl;
										break;
									}
					}
					
					t += timeTau;
				}
										
				

				//cout<<"FIRED REACTION: "<<idOfReac<<" and time is "<<t<<endl;
			}
			
			//auto eIter = chrono::high_resolution_clock::now() - aIter;
			//cout << "Iter: " << chrono::duration_cast<chrono::nanoseconds>(eIter).count() << " nanoseconds"<< endl; 
			for(int k=0; k<M; k++)
				simulation[k] += x[k];
			
			auto eSim = chrono::high_resolution_clock::now() - aSim;
			//cout << "total time for Rand: " << eRand << " nanoseconds"<< endl;
			//cout << "total time for 1 Sim: " << chrono::duration_cast<chrono::nanoseconds>(eSim).count() << " nanoseconds"<< endl; 
			//cout << "time of 1 sim without rand: " << chrono::duration_cast<chrono::nanoseconds>(eSim).count() - eRand << " nanoseconds"<< endl; 
			//cout<<"used rand numbers="<<count<<endl;
		}

		//cout<<"Sending to collector"<<endl;
		//for(int i=0; i<M; i++)
		//	cout<<simulation[i]<<" ";
		//cout<<endl;
	
		return (new vector<int>(simulation));
		
	};

	void eosnotify(ssize_t) 
	{
			ff_send_out(EOS);	
	};

	void init()
	{		
		// reaction rates	
		rate[0] = 1;		 	withDelay[0] = 2; //"NCD";
		rate[1] = 0.01;		 	withDelay[1] = 3; //"ND";
		rate[2] = 5;		 	withDelay[2] = 3; //"ND";
		rate[3] = 0.01;		 	withDelay[3] = 1; //"CD";
		
		
		// S1 --r0--> S1 + S2	
		reactant[0][0] = -1;
		product[0][0] = +1;
		product[0][1] = +1;
		
		
		// S1 + S2 --r1--> S3
		reactant[1][0] = -1;
		reactant[1][1] = -1;
		product[1][2] = +1;
		
		// S3 --r2--> S1 + S2
		reactant[2][2] = -1;
		product[2][0] = +1;
		product[2][1] = +1;
				
		// S2 --r3--> _
		reactant[3][1] = -1;
		
	};	
	
	void generateRandomArrays(double seed, int __npf)
	{
		ff::ParallelFor pf(__npf);
		//cout<<"here"<<endl;
		
		pf.parallel_for(0L, arraysize, 1, 0, [&](const long i)
		{
			myrandQueue[i] = randDouble(seed);
			//cout<<"myrandQueue="<<myrandQueue[i]<<endl;
		});
	}
	
	double randDouble(double seed) 
	{
		//thread_local static std::random_device rd;
		thread_local static std::mt19937 rng(seed);
		thread_local std::uniform_real_distribution<double> urd;
		return urd(rng, decltype(urd)::param_type{0,1});
	}

	bool Check_reaction(int idOfReac)
	{
		for(int j=0; j<M; j++)
		{
			if(reactant[idOfReac][j]<0 && x[j]<reactant[idOfReac][j]*(-1))
			{
				//cout<<"there are enough reactant to fire this reaction"<<endl;
				return false;
			}
		}	
		return true;
	}	
	
	double propensity_sum()
	{
		p[0] = rate[0]*x[0];		    		// S1 --r0--> S1 + S2
		p[1] = rate[1]*x[0]*x[1];				// S1 + S2 --r1--> S3
		p[2] = rate[2]*x[2];		    		// S3 --r2--> S1 + S2
		p[3] = rate[3]*x[1];		    		// S2 --r3--> _		    
	
		double p_sum=0.0;
		
		for(int i=0; i<N; i++)
			p_sum += p[i];
		
		return p_sum;
	};

	
};

// Third Stage. Collector of the Farm
struct Collector: ff_node_t< vector<int>, int > 
{		
	int num_simulation;
	int simulation[M];
	
	Collector(int num_simulation):num_simulation(num_simulation){}

	int svc_init()
	{
		for(int i=0; i<M; i++)
			simulation[i]=0;
		return 0;
	}		
	
	int* svc(vector<int> *x) 
	{				
		vector<int> &v = *x;
		for(int i=0; i<M; i++)
		{
			simulation[i] += v[i];
		}
		//cout<<endl;
		return GO_ON;		
    };
	
	void svc_end()
	{	
		for(int i=0; i<M; i++)
		{
			simulation[i] = simulation[i] / num_simulation;
		}
		
		for(int i=0; i<M; i++)
			cout<<simulation[i]<<" ";
		cout<<endl;		
	};	

};

int main(int argc, char * argv[])
{
	 if (argc<4) 
	 {
        std::cerr << "use: " << argv[0]  << " number1 number2 number3 as num_worker num_simulation num_par_for\n";
        return -1;
    }
	
	int num_worker 		= std::stol(argv[1]);	// The number of workers
	int num_simulation	= std::stol(argv[2]);	// The number of simulations
	int num_par_for		= std::stol(argv[3]);	// The number of workers for parallel for
	
//************** Creating nodes of FastFlow *******************
	
	Emitter emitter(num_worker);	
	Collector collect(num_simulation);
	
	std::vector<std::unique_ptr<ff_node> > W; 
	for(int i=0;i<num_worker;++i) 
		W.push_back(make_unique<Simulation>(num_worker, num_simulation, num_par_for));
											
	ff_Farm<int> farm(std::move(W), emitter, collect); 				
	ff_Pipe<> pipe(farm);
	
//*************************************************************
	
	
//************************** Fast Flow ************************	
	auto startReac = chrono::high_resolution_clock::now();
	ffTime(START_TIME);
	
	if (pipe.run_and_wait_end()<0) 
	{
		error("running pipe");
		return -1;
	}	
	
	auto elapsedReac = chrono::high_resolution_clock::now() - startReac;
	ffTime(STOP_TIME);
	
	
	cout<<"Time	ffTime(GET_TIME): "<< ffTime(GET_TIME) <<endl;		
	cout << "Time	chrono: " << chrono::duration_cast<chrono::microseconds>(elapsedReac).count() << " microseconds"<< endl; 
	//std::cout << "Time	pipe.ffTime(): " << farm.ffTime() << "\n";	
//*****************************************************************
	
	
	
	
	
//pipe.ffStats(std::cout);
	return 0;
}

