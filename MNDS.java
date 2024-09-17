import java.util.Arrays;
import java.util.BitSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class MNDS {
    public static void main(String[] args) {
	/* int n = 4; 
	int d = 3;
	Point[] population = new Point[n];
		
        double[] objectives1 = new double[3];
        objectives1[0] = 2; objectives1[1] = 3; objectives1[2] = 1;
        population[0] = new Point(0, n, objectives1);
        
        double[] objectives2 = new double[3];
        objectives2[0] = 1; objectives2[1] = 4; objectives2[2] = 1;
        population[1] = new Point(1, n, objectives2);
        
        double[] objectives3 = new double[3];
        objectives3[0] = 3; objectives3[1] = 2; objectives3[2] = 1;
        population[2] = new Point(2, n, objectives3);
        
        double[] objectives4 = new double[3];
        objectives4[0] = 4; objectives4[1] = 1; objectives4[2] = 1;
        population[3] = new Point(3, n, objectives4); */





	/* int n = 4; 
	int d = 3;
	Point[] population = new Point[n];
		
        double[] objectives1 = new double[3];
        objectives1[0] = 11; objectives1[1] = 12; objectives1[2] = 13;
        population[0] = new Point(0, n, objectives1);
        
        double[] objectives2 = new double[3];
        objectives2[0] = 13; objectives2[1] = 12; objectives2[2] = 11;
        population[1] = new Point(1, n, objectives2);
        
        double[] objectives3 = new double[3];
        objectives3[0] = 4; objectives3[1] = 5; objectives3[2] = 6;
        population[2] = new Point(2, n, objectives3);
        
        double[] objectives4 = new double[3];
        objectives4[0] = 6; objectives4[1] = 5; objectives4[2] = 4;
        population[3] = new Point(3, n, objectives4);*/
        
		
		     
		
		
        /* int n = 6; 
	int d = 4;
        Point[] population = new Point[n];
        
        double[] objectives1 = new double[4];
        objectives1[0] = 1; objectives1[1] = 6; objectives1[2] = 1; objectives1[2] = 1;
        population[0] = new Point(0, n, objectives1);
        
        double[] objectives2 = new double[4];
        objectives2[0] = 2; objectives2[1] = 1; objectives2[2] = 2; objectives2[2] = 3;
        population[1] = new Point(1, n, objectives2);
        
        double[] objectives3 = new double[4];
        objectives3[0] = 5; objectives3[1] = 2; objectives3[2] = 4; objectives3[2] = 5;
        population[2] = new Point(2, n, objectives3);
        
        double[] objectives4 = new double[4];
        objectives4[0] = 4; objectives4[1] = 4; objectives4[2] = 6; objectives4[2] = 4;
        population[3] = new Point(3, n, objectives4);
        
        double[] objectives5 = new double[4];
        objectives5[0] = 3; objectives5[1] = 5; objectives5[2] = 2; objectives5[2] = 3;
        population[4] = new Point(4, n, objectives5);
        
        double[] objectives6 = new double[4];
        objectives6[0] = 6; objectives6[1] = 3; objectives6[2] = 6; objectives6[2] = 8;
        population[5] = new Point(5, n, objectives6); */

        
        Helper helper = new Helper();
        int n = 10;
        int d = 4;
        long seed = 44;
        Point[] population = new Point[n];
        population = helper.geneartedata(n, d, seed); // Generate a random population

        helper.printPopulation(population);
        
        helper.sortMNDS(population); // Sort the population using Merge Non-dominated Sort
    }
}


class Point {
    private int id;
    private double[] objectives;
    private BitSet dominanceSet;
    private int index = -1;

    public Point(int id, int populationSize) {
        this.id = id;
        this.dominanceSet = new BitSet(populationSize);
    }

    public Point(int id, int populationSize, double[] objectives) {
        this.id = id;
        this.dominanceSet = new BitSet(populationSize);
        this.objectives = objectives;
    }
    
    public int getId() {
        return this.id;
    }

    public double[] getObjectives() {
        return this.objectives;
    }

    public BitSet getDominanceSet() {
        return this.dominanceSet;
    }

    public int getIndex() {
        return index;
    }

    public void setId(int id) {
        this.id = id;
    }

    public void setObjectives(double[] objectives) {
        this.objectives = objectives;
    }

    public void setDominanceSet(BitSet dominanceSet) {
        this.dominanceSet = dominanceSet;
    }

    public void set(int bitIndex) {
        this.dominanceSet.set(bitIndex);
    }

    public void setIndex(int index) {
        this.index = index;
    }

    /*  1: objectives dominates point
       -1: point dominates objectives
        0: objectives and point is non-dominated */
    public int dominanceRelationship(Point p) {
        return dominanceComparison(p.objectives, this.objectives, this.objectives.length);
    }
    
    public int dominanceComparisonHelper(double[] a, double[] b, int dim, int maxDim, int returnOnEnd) {
        while (++dim < maxDim) {
            if (a[dim] < b[dim]) {
                return 0;
            }
        }
        return returnOnEnd;
    }

    public int dominanceComparison(double[] a, double[] b, int dim) {
        for (int i = 0; i < dim; ++i) {
            if (a[i] < b[i]) {
                return dominanceComparisonHelper(b, a, i, dim, -1);
            }
            if (a[i] > b[i]) {
                return dominanceComparisonHelper(a, b, i, dim, 1);
            }
        }
        return 0;
    }
    
    /* Used for pre-sorting as in ENS
     * 1: First point is having small value for a objctive function 
      -1: Second point is having small value for a objctive function 
       0: Same */
    public int isSmall(Point p) {
        double[] a = p.objectives;
        int noObjectives = a.length;
        for(int i = 0; i < noObjectives; i++) {
            if(this.objectives[i] < a[i]) {
                return 1;
            } else if (this.objectives[i] > a[i]) {
                return -1;
            }
        }
        return 0;
    }
    
    public int isSmall(Point p, int[] Q0tOrder, int m) {
        if(this.objectives[m] < p.objectives[m]) {
            return 1;
        } else if (this.objectives[m] > p.objectives[m]) {
            return -1;
        } else {
            if(Q0tOrder[this.id] < Q0tOrder[p.id]) {
                return 1;
            } else {
                return -1;
            }
        }
    }
    
    public boolean isEqual(Point p) {
        double[] a = p.objectives;
        int noObjectives = a.length;
        for(int i = 0; i < noObjectives; i++) {
            if(this.objectives[i] != a[i]) {
                return false;
            } 
        }
        return true;
    }

    @Override
    public String toString() {
        return "Point{" + "id=" + id + ", objectives=" + Arrays.toString(objectives) + ", dominanceSet=" + dominanceSet + ", index=" + index + '}';
    }
}


class Helper {
    public Point[] geneartedata(int n, int d, long seed) {
        Random random = new Random(seed);
        Point[] population = new Point[n];
        for (int i = 0; i < n; ++i) {
            population[i] = new Point(i, n);
            double[] point = new double[d];
            for (int j = 0; j < d; ++j) {
                point[j] = random.nextDouble();
            }
            population[i].setObjectives(point);
        }
        return population;
    }
    
    public void sortMNDS(Point population[]) {
        int n = population.length; // Number of points in the population
        int noObj = population[0].getObjectives().length; // Size of the objective vector associated with the point
  
        BitSet ods = new BitSet(n);
        int odsMin = n; // The index of the minimum set bit in the ods bitset; It can be obtained using 'ods.nextSetBit(0)' function
        int odsMax = -1; // The index of the maximum set bit in the ods bitset
        
        HeapSort hs = new HeapSort();
        int[] Q0 = new int[n]; // Stores the ids of the sorted points based on the first objective 
        int[] Q1 = new int[n]; // Stores the ids of the sorted points based on the second objective 
        int[] Qm = new int[n]; // Stores the ids of the sorted points based on the remaining objectives
        int[] ranks = new int[n]; // Array to store the rank of the points
        int PIndex;
        for(int i = 0; i < n; i++) {
            Q0[i] = population[i].getId();
            Q1[i] = population[i].getId();
            Qm[i] = population[i].getId();
            ranks[i] = 1; // Initialize the rank of all the points to 1
        }

        Q0 = hs.sort1(population, Q0); // Sort the points based on the first objective
        for(int i = 0; i < n; i++) {
            population[Q0[i]].setIndex(i);
        }
        System.out.println("\nQ0 = " + Arrays.toString(Q0));
        printPopulation(population);

        Q1 = hs.sortm(population, Q1, Q0, 1); // Sort the points based on the second objective
        System.out.println("\nQ1 = " + Arrays.toString(Q1));
        printPopulation(population);
        
        /* Obtain the dominance set based on the first two objectives */
        boolean hasDominance = false; // A boolean variable to chcek whether any point is dominated based of the first two objectives
        for(int i = 0; i < n; i++) {
            PIndex = population[Q1[i]].getIndex();
            if(PIndex < odsMin) { // Subset is empty
                odsMin = PIndex;
            } else {
                //population[Q2[i]].getDominanceSet().and(ods);
                /* Subset method */
                int top = Math.min(PIndex-1, odsMax);
//                for (int j = odsMin; j <= top && j >= 0; j = ods.nextSetBit(j+1)) {
                for (int j = odsMin; j >= 0; j = ods.nextSetBit(j+1)) {
                    population[Q1[i]].set(j);
                }
                hasDominance = true;
            }
            ods.set(PIndex);
            if(PIndex > odsMax) {
                odsMax = PIndex;
            }
        }
        
        if(hasDominance) { // If one of the point is dominated considering the dominance set based on the first two objectives 
            int m = 2;
            while(m < noObj && hasDominance) {
                System.out.println("\nm = " + m);
                Qm = hs.sortm(population, Qm, Q0, m); // Sort the points based on the m-th objective
                System.out.println("Qm = " + Arrays.toString(Qm));
                hasDominance = false;
                ods = new BitSet(n);
                for(int i = 0; i < n; i++) {
                    PIndex = population[Qm[i]].getIndex();
                    population[Qm[i]].getDominanceSet().and(ods); 
                    ods.set(PIndex);
                    if(!population[Qm[i]].getDominanceSet().isEmpty()) {
                        hasDominance = true; // Point Qm[i] is dominated so hasDominance is set to True
                    }
                }
                printPopulation(population);
                m = m + 1;
            }
            
            int[] QIndex = new int[n];
            for(int i = 0; i < n; i++) {
                QIndex[i] = population[Qm[i]].getIndex();
            }

            ranks[Qm[0]] = 1;
            for(int i = 1; i < n; i++) {
                int myRank = 1;
                for (int j = population[Qm[i]].getDominanceSet().nextSetBit(0); j >= 0; 
                    j = population[Qm[i]].getDominanceSet().nextSetBit(j+1)) {
                    int thatRank = ranks[Qm[QIndex[j]]];
                    if(thatRank >= myRank) {
                        myRank = thatRank + 1;
                    }
                }
                ranks[Qm[i]] = myRank;
            }
            System.out.println("ranks = " + Arrays.toString(ranks));
        } else { /* If no point is dominated considering the dominance set based on the first two objectives. 
            In this case, all the points belong to first fornt and thus having rank 1 */
            System.out.println("ranks = " + Arrays.toString(ranks));
        }
    }
    
    public void printPopulation(Point[] population) {
        for (Point pop : population) {
            System.out.println(pop);
        }
    }
    
    public void sortDeductiveSort(Point population[]) {
        System.out.println("************* DS ************* ");
        double noDC = 0; // No of dominance comparisons
        int timecomp = 0;
        double noCheck = 0; // No of dominance comparisons
        
        int x = 0;
        int f = 0;
        ArrayList<ArrayList<Integer>> setNonDominatedFront = new ArrayList();
        boolean[] isSorted = new boolean[population.length];
        
        while(f < population.length) { 
            boolean[] D = new boolean[population.length];
            ArrayList<Integer> ndf = new ArrayList();
            for(int i = 0; i < population.length; i++) {
                int count = 0;
                if(!D[i] && !isSorted[i]) {
                    for(int j = i+1; j < population.length; j++) {
                        count++;
                        if(!D[j] && !isSorted[j]) {
                            int isDominates = population[i].dominanceRelationship(population[j]);
                            noDC++;
                            if(isDominates == 1) { // i dominates j
                                D[j] = true;
                            } else if (isDominates == -1) { // j dominates i
                                D[i] = true;
                                break;
                            }
                        }
                    }
                    timecomp = timecomp + count;
                    if(!D[i]) {
                        ndf.add(population[i].getId());
                        isSorted[i] = true;
                        f++;
                    }
                } else {

                }
            }
            setNonDominatedFront.add(ndf);
        }
        
        for(ArrayList<Integer> element : setNonDominatedFront) {
            Collections.sort(element);
        }
        System.out.println("The sorted fronts....");
        for(ArrayList<Integer> element : setNonDominatedFront) {
            System.out.println(element);
        }
    }
}


class HeapSort {
    void heapifyFirstObj(int heapSize, int i, Point population[], int Q[]) {
        int largest = i;  // Initialize largest as root
        int l = 2*i + 1;  // left = 2*i + 1
        int r = 2*i + 2;  // right = 2*i + 2
         
        if (l < heapSize && population[Q[l]].isSmall(population[Q[largest]]) == -1) {
                largest = l;
        }

        if (r < heapSize && population[Q[r]].isSmall(population[Q[largest]]) == -1) {
                largest = r;
        }

        if (largest != i) {
            int swap = Q[i];
            Q[i] = Q[largest];
            Q[largest] = swap;
            heapifyFirstObj(heapSize, largest, population, Q);
        }
    }
    
    void heapifyOtherObj(int Q[], int heapSize, int i, Point population[], int[] Q0tOrder, int m) {
        int largest = i;  // Initialize largest as root
        int l = 2*i + 1;  // left = 2*i + 1
        int r = 2*i + 2;  // right = 2*i + 2

        if (l < heapSize && population[Q[l]].isSmall(population[Q[largest]],Q0tOrder, m) == -1) {
            largest = l;
        }

        if (r < heapSize && population[Q[r]].isSmall(population[Q[largest]],Q0tOrder, m) == -1) {
            largest = r;
        }

        if (largest != i) {
            int swap = Q[i];
            Q[i] = Q[largest];
            Q[largest] = swap;
            heapifyOtherObj(Q, heapSize, largest, population, Q0tOrder, m);
        }
    }
    
    public void print(Point population[]) {
        for (Point point : population) {
            System.out.println(point);
        }
    }
    
    /*-- Sort based on first objective --*/
    public int[] sort1(Point population[], int Q[]) {
        int n = Q.length;
        for (int i = n / 2 - 1; i >= 0; i--) {
            heapifyFirstObj(n, i, population, Q);
        }
        for (int i=n-1; i>=0; i--) {
            int temp = Q[0];
            Q[0] = Q[i];
            Q[i] = temp;
            heapifyFirstObj(i, 0, population, Q);
        }
        return Q;
    }
    
    /*-- Sort based on other objective --*/
    public int[] sortm(Point population[], int Q[], int Q0tOrder[], int m) {
        int n = Q.length;
        for (int i = n / 2 - 1; i >= 0; i--) {
            heapifyOtherObj(Q, n, i, population, Q0tOrder, m);
        }
        for (int i=n-1; i>=0; i--) {
            int temp = Q[0];
            Q[0] = Q[i];
            Q[i] = temp;
            heapifyOtherObj(Q, i, 0, population, Q0tOrder, m);
        }
        return Q;
    }
}
