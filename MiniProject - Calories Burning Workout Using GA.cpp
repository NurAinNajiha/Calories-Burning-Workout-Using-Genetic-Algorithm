/* EVO Mini Project
Title : Calories-Burning Workout Using  Genetic Algorithm
TeamMates : Azfar Rahman Bin Fazul Rahman (B031910467) - Programmer,algorithm designer, system analyst
			Nur’Ain Najiha binti Zakaria (B031910466) -  Documentation,algorithm designer
			Nur Izzati binti Shafie      (B031910476) -  Documentation,algorithm designer
			Megala d/o Sontulom          (B031910172) -  Documentation,algorithm designer
*/

//libaries
#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;


//Declaration variable and representation of a chromosome
const int CAPACITIES_CALORIES = 1000;
const int GENE = 15;
const int POP_SIZE = 30;// increased from 10 to 30  
const int CALORIES[GENE] = {/*Skipping 15mins*/300,/*jogging 10mins*/150,/*jumping jacks 1min*/16,/* Squats 25mins*/ 97,/* Plank 1min*/ 5,/* push up 1min*/ 7,/* pull up 1min*/ 7,/*Walking Lunges 10min*/ 80,/* Swimming 1hour*/ 400,/* cycling 30mins*/ 298,/* weightlifting 30mins*/ 126,/* Aerobic 30mins*/ 295,/* Abdominal Crunches 1mins*/ 3,/*stretching  30mins*/ 70,/* yoga 30mins*/ 180 };
const int TIME[GENE] = { 15,10,1,25,1,1,1,10,60,30,30,30,1,30,30 };//time taken for every exercise 
const float CROSS_PROB = 0.2; //Probability for the uniform crossover to occur
const float MUT_PROB = 0.1;//Probability for the mutation
const int GENERATION = 100;//max generation
const int K = 14;// number of players in tournement selection
const float coinflipProbability = 0.5; // probability for each gene while in crossover


//declaration of chromosome, fitness,parent,children,newChromosome ,totalrun
//array in the program
int chromosome[POP_SIZE][GENE];
double fitness[POP_SIZE];
int parents[2][GENE];
int children[2][GENE];
int newChromo[POP_SIZE][GENE];
double Totalruning[POP_SIZE];

//dclaration of best chromosome,best fitness,average fitness
int bestChromo[GENE];
double bestFit = 99999999999.9;
double avgFit = 0.0;

//file
ofstream BF, AF, BC;

//initializing the popultion
void initializePopulation() {
	
	srand(time(NULL));
	for (int c = 0; c < POP_SIZE; c++) {
		for (int i = 0; i < GENE; i++) {
			
			chromosome[c][i] = rand() % 2;
		}
	}
}

void printChromo() {
	for (int c = 0; c < POP_SIZE; c++) {


		cout << "\tC" << c << ": \t";
		for (int i = 0; i < GENE; i++) {
			cout << chromosome[c][i] << " ";
		}
		cout << endl;
	}
}

//Evaluating the chromosomes.
void evaluateChromo() {
	cout << "\t Chromosome\tDifference\tTotal Time (Minutes)\tFitness Value (FV)" << endl;
	int accumulatedCalories, accumulatedTime;
	for (int c = 0; c < POP_SIZE; c++) {
		accumulatedCalories = 0;
		accumulatedTime = 0;
		for (int i = 0; i < GENE; i++) {
			if (chromosome[c][i] == 1) {
				accumulatedCalories = accumulatedCalories + CALORIES[i];
				accumulatedTime = accumulatedTime + TIME[i];
			}
		}
		fitness[c] = abs(CAPACITIES_CALORIES - (accumulatedCalories + accumulatedTime)) / (float)CAPACITIES_CALORIES;
		cout << "\t C" <<c<< "\t\t" << abs(CAPACITIES_CALORIES - accumulatedCalories) << "\t\t " << accumulatedTime << "\t\t\t" << fitness[c] << endl;
	}
}



//parent : Tournament Selection
void TourSelection() {

	int player[K] = {};
	int ParentsIndex[2];
	double PlayerFit[K] = {};
	int IndexTaken[K];
	int mostFitIndex;
	int SecondmostFitIndex;
	double mostFitSofar = 9999999;
	double secondMostFit = 9999;
    
    srand(time(NULL));
	//players
StartAgain:
	do{
	
	for (int p = 0; p < K; p++) {
		player[p] = rand() % (POP_SIZE);
		for (int e = 0; e < p; e++) {
			if (player[p] == player[e])
				goto StartAgain;
		}
		IndexTaken[p] = player[p];

	}

	//this is for the FIRST ROUND 
	for (int i = 0; i < K; i++) {
		PlayerFit[i] = fitness[IndexTaken[i]];
		if (PlayerFit[i] < mostFitSofar) {
			mostFitSofar = PlayerFit[i];
			mostFitIndex = IndexTaken[i];
		}
	}


	//this is for the SECOND ROUND
	for (int j = 0; j < K; j++) {
		if (IndexTaken[j] == mostFitIndex)
			continue;
		PlayerFit[j] = fitness[IndexTaken[j]];
		if (PlayerFit[j] < secondMostFit) {
			secondMostFit = PlayerFit[j];
			SecondmostFitIndex = IndexTaken[j];
		}
	}
}
while(IndexTaken[0]==IndexTaken[1]);


	//the output should be
	cout << "\tPlayers\t\tFV" << endl;
	for (int d = 0; d < K; d++) {
		cout << "\tC" << IndexTaken[d] << "\t\t" << PlayerFit[d] << endl;
	}
	cout << "\n\tMost Fit\tFV" << endl;
	cout << "\tC" << mostFitIndex << "\t\t" << mostFitSofar << "\n\tC" << SecondmostFitIndex << "\t\t" << secondMostFit << endl;


	//then insert to parents chromo
	cout << "\n\tParents\t\tChromosome" << endl;
	for (int p = 0; p < 2; p++) {
		//parent 1
		if (p == 0) {
			cout << "\tC" << mostFitIndex << "\t\t";
			for (int g = 0; g < GENE; g++) {
				parents[p][g] = chromosome[mostFitIndex][g];
				cout << chromosome[mostFitIndex][g] << " ";
			}
			cout << endl;
		}

		//parent 2
		if (p == 1) {
			cout << "\tC" << SecondmostFitIndex << "\t\t";
			for (int g = 0; g < GENE; g++) {
				parents[p][g] = chromosome[SecondmostFitIndex][g];
				cout << chromosome[SecondmostFitIndex][g] << " ";

			}
			cout << endl;
		}
	}
}



//Uniform Crossover
void crossover() {

	//in this crossover we have flip a coin in each gene of a parent.

	float prob;
	int CrossOverPoint;
	float CrossoverProb;

	//checking if crossover will happen

	CrossoverProb = ((rand()) % 10) + 1 / 10.0;

	if (CrossoverProb > CROSS_PROB) {
		////Crossover did not happen
		cout << "\t Crossover not happened.Parents will be copied to children" << endl;
		//copying from parent to children 
		for (int p = 0; p < 2; p++) {
			for (int g = 0;g < GENE; g++) {
				children[p][g] = parents[p][g];
			}
			cout << endl;
		}
		cout << "\t******************************************" << endl;
		//children Output
		for (int c = 0; c < 2; c++) {
			cout << "\tChildren" << c + 1 << "\t";
			for (int g = 0; g < GENE; g++) {
				cout << children[c][g] << " ";

			}
			cout << endl;
		}
	}
	else {
		//crossover happens.
		//copying parents to children.
		for (int p = 0; p < 2; p++) {
			for (int g = 0; g < GENE; g++) {
				children[p][g] = parents[p][g];
			}
		}
		//flipping coin in each gene, it wil be exchange based on #coinflipProbability.
		for (int chromoCount = 0; chromoCount < GENE;chromoCount++) {
			prob = ((rand() % 10) + 1) / 10.0;
			if (prob <= coinflipProbability) {
				CrossOverPoint = chromoCount;
				children[0][chromoCount] = parents[1][chromoCount];
				children[1][chromoCount] = parents[0][chromoCount];
				cout << "\tExchange gene " << chromoCount + 1 << endl;
			}
		}
		cout << "\tDONE CROSSOVER!" << endl << endl;

		//parents original output
		for (int c = 0; c < 2; c++) {
			cout << "\tParent" << c + 1 << ":\t";
			for (int i = 0; i < GENE; i++) {
				cout << parents[c][i] << " ";
			}
			cout << endl;
		}
		cout << "\t******************************************" << endl;
		//children Output
		for (int c = 0; c < 2; c++) {
			cout << "\tChildren" << c + 1 << ":\t";
			for (int g = 0; g < GENE; g++) {
				cout << children[c][g] << " ";
			}
			cout << endl;
		}
	}
}

// mutation process
void mutation() {
	float prob;
	int mut_point;

	for (int c = 0; c < 2; c++) {
		prob = (rand() % 11) / 10.0;
		if (prob < MUT_PROB) {
			mut_point = rand() % GENE;
			cout << "\tMutation at gene " << mut_point + 1 << " for child " << c + 1 << endl;
			if (children[c][mut_point] == 0)
				children[c][mut_point] = 1;
			else
				children[c][mut_point] = 0;
		}
		else
			cout << "\tMutation Never Occured for children : " << c + 1 << endl;
	}
	cout << "\t********************************************";
	for (int c = 0; c < 2; c++) {
		cout << "\n\tChildren " << c + 1 << " after mutation:\t";
		for (int g = 0; g < GENE; g++) {
			cout << children[c][g] << " ";
		}
	}
}


// Survival slection based on FItness-Based Selection
void SurvivalSelection() {
	//take from population - actually take from mutation.
	// take from K players from population
	// Choose 2 winners - 2 winners have the highest fitness (bad)
	// Replace the 2 winners with the children AFTER mutation.
	
	
    //In this fitness based selection, the children will replace
    //the least fit individuals in the population. In this case,
    //we will select the least fit individuals using tournament selection.
    //1. Select Least Fit Individuals using K-way Tournament Selection.

	int player[K] = {};
	int Replacedindex[2] = {};
	double playerFit[K] = {};
	int Indextaken[K] = {};
	int leastFitIndex;
	int secondLeastFitIndex;
	double leastFitSoFar = 0;
	double secondLeastFit = 0;


	// Gettting players for the tournament.
	for (int p = 0; p < K; p++) {
		player[p] = rand() % (POP_SIZE);
		Indextaken[p] = player[p];
	}

	// process of calculating fitness
	//FIRST ROUND
	for (int i = 0; i < K; i++) {
		playerFit[i] = fitness[Indextaken[i]];
		if (playerFit[i] > leastFitSoFar) {
			leastFitSoFar = playerFit[i];
			leastFitIndex = Indextaken[i];
		}
	}
	// Second Round
	for (int j = 0; j < K; j++) {
		if (Indextaken[j] == leastFitIndex)
			continue;
		playerFit[j] = fitness[Indextaken[j]];
		if (playerFit[j] > secondLeastFit) {
			secondLeastFit = playerFit[j];
			secondLeastFitIndex = Indextaken[j];
		}
	}

	// Output
	cout << "\tPlayers\t\tFV" << endl;
	for (int d = 0; d < K; d++) {
		cout << "\tC" << Indextaken[d] << "\t\t" << playerFit[d] << endl;
	}
	cout << "\n\tLeast Fit\tFV" << endl;
	cout << "\tC" << leastFitIndex << "\t\t" << leastFitSoFar << "\n\tC" << secondLeastFitIndex << "\t\t" << secondLeastFit << endl;

	// Prepare Buffer Chromosomes
	for (int c = 0; c < POP_SIZE; c++) {
		for (int g = 0; g < GENE; g++) {
			newChromo[c][g] = chromosome[c][g];
		}
	}

	// Copy children to replace the parents.
	for (int c = 0; c < 2; c++) {
		if (c == 0) {
			for (int g = 0; g < GENE; g++) {
				newChromo[leastFitIndex][g] = children[c][g];
			}
		}
		else {
			for (int g = 0; g < GENE; g++) {
				newChromo[secondLeastFitIndex][g] = children[c][g];
			}
		}
	}
	cout << endl;

	// Copy the chromosomes from buffer to population
	for (int c = 0; c < POP_SIZE; c++) {
		for (int g = 0; g < GENE; g++) {
			chromosome[c][g] = newChromo[c][g];
		}
	}

}


void recordBestFitness() {
   
	/*
	int bestChromosomeIndex;
    for (int c = 0; c < POP_SIZE; c++) {
        if (bestFit > fitness[c]) {
            bestFit = fitness[c];
            bestChromosomeIndex = c;
            for (int g = 0; g < GENE; g++) {
                bestChromo[g] = chromosome[c][g];
            }
        }
    }

    cout << "\n\tBest Fitness = " << bestFit;
    cout << "\n\tBest Chromosome: ";
    for (int g = 0; g < GENE; g++) {
        cout << bestChromo[g] << " ";
    }
    cout << endl;

    // Write to text file
    BF << bestFit << endl;
    for (int g = 0; g < GENE; g++) {
        BC << bestChromo[g] << " ";
    }
    BC<< endl;
	*/
	for (int i = 0;i < POP_SIZE;i++)
	{
		if (bestFit > fitness[i])
		{
			bestFit = fitness[i];

			for (int g = 0;g < GENE;g++)
			{
				bestChromo[g] = chromosome[i][g];
			}
				
		}
	}

	cout << "\nBEST FITNESS : " << bestFit << endl;
	cout << "\nBEST CHROMOSOME : ";

	for (int i = 0;i < GENE;i++)
	{
		cout << bestChromo[i];
	}

	cout << endl;

	BF << bestFit << endl;

	for (int i = 0;i < GENE;i++)
	{
		BC << bestChromo[i] << " ";
	}

	BC << endl;
}

void calculateAvgFitness()
{
	/*
	double sum = 0;

	for (int i = 0;i < POP_SIZE;i++)
	{
		sum = sum + fitness[i];
	}
	avgFit = sum / POP_SIZE;

	cout << "\nAVERAGE FITNESS : " << avgFit << endl;
	AF << avgFit << endl;

	*/
	double sum = 0;

	for (int i = 0;i < POP_SIZE;i++)
	{
		sum = sum + fitness[i];
	}
	avgFit = sum / POP_SIZE;

	cout << "\nAVERAGE FITNESS : " << avgFit << endl;
	AF << avgFit << endl;
}

int main() {
	//initialization
srand(time(0));
	BC.open("BestChromo.txt");
	BF.open("BestFitness.txt");
	AF.open("AvgFitness.txt");
	cout << "\tLets Start The Genatic Algorithm." << endl << endl;
	cout << "\tInitial Population" << endl;
	initializePopulation();
	printChromo();

	for (int g = 0; g < GENERATION; g++) {
		cout << "\t******** GENERATION" << g + 1 << " ************" << endl;
		//Evaluation
		cout << endl;
		cout << "\t1. Evaluation of chromosomes" << endl;
		evaluateChromo();
		recordBestFitness();
		calculateAvgFitness();

		//parent selection 
		cout << "\n\t2. PARENT SELECTION : " << K << "-way Tournament" << endl;
		TourSelection();


		//uniform crossover
		cout << endl;
		cout << "\n\t3. Crossover: UNIFORM CROSSOVER PROCESS" << endl;
		crossover();
		cout << endl;


		//Mutation
		cout << endl;
		cout << "\n\t3. MUTATION PROCESS" << endl;
		mutation();
		cout << endl;


		// Survival Selection
		cout << endl;
		cout << "\n\t5. Survival Selection: " << K << "-way Tournament Selection and Fitness Based Selection" << endl;
		SurvivalSelection();

	}
	//bestFitnessFile
	BF.close();
	//avgFitnessFile
	AF.close();
	//bestChromosomeFile
	BC.close();
	cout << endl;
	
}


