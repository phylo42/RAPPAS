/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos.sorting;

import core.ProbabilisticWord;
import java.util.List;

/**
 *
 * @author ben
 */
public class HeapSort {
    
    public static void heapSort(int[] a){
	int count = a.length;
 
	//first place a in max-heap order
	heapify(a, count);
 
	int end = count - 1;
	while(end > 0){
		//swap the root(maximum value) of the heap with the
		//last element of the heap
		int tmp = a[end];
		a[end] = a[0];
		a[0] = tmp;
		//put the heap back in max-heap order
		siftDown(a, 0, end - 1);
		//decrement the size of the heap so that the previous
		//max value will stay in its proper place
		end--;
	}
    }

    private static void heapify(int[] a, int count){
            //start is assigned the index in a of the last parent node
            int start = (count - 2) / 2; //binary heap

            while(start >= 0){
                    //sift down the node at index start to the proper place
                    //such that all nodes below the start index are in heap
                    //order
                    siftDown(a, start, count - 1);
                    start--;
            }
            //after sifting down the root all nodes/elements are in heap order
    }

    private static void siftDown(int[] a, int start, int end){
            //end represents the limit of how far down the heap to sift
            int root = start;

            while((root * 2 + 1) <= end){      //While the root has at least one child
                    int child = root * 2 + 1;           //root*2+1 points to the left child
                    //if the child has a sibling and the child's value is less than its sibling's...
                    if(child + 1 <= end && a[child] < a[child + 1])
                            child = child + 1;           //... then point to the right child instead
                    if(a[root] < a[child]){     //out of max-heap order
                            int tmp = a[root];
                            a[root] = a[child];
                            a[child] = tmp;
                            root = child;                //repeat to continue sifting down the child now
                    }else
                            return;
            }
    }
    
    //-------------------------------------------------------------------
    
    public static void heapSort(double[] a){
	int count = a.length;
 
	//first place a in max-heap order
	heapify(a, count);
 
	int end = count - 1;
	while(end > 0){
		//swap the root(maximum value) of the heap with the
		//last element of the heap
		double tmp = a[end];
		a[end] = a[0];
		a[0] = tmp;
		//put the heap back in max-heap order
		siftDown(a, 0, end - 1);
		//decrement the size of the heap so that the previous
		//max value will stay in its proper place
		end--;
	}
    }

    private static void heapify(double[] a, int count){
            //start is assigned the index in a of the last parent node
            int start = (count - 2) / 2; //binary heap

            while(start >= 0){
                    //sift down the node at index start to the proper place
                    //such that all nodes below the start index are in heap
                    //order
                    siftDown(a, start, count - 1);
                    start--;
            }
            //after sifting down the root all nodes/elements are in heap order
    }

    private static void siftDown(double[] a, int start, int end){
            //end represents the limit of how far down the heap to sift
            int root = start;

            while((root * 2 + 1) <= end){      //While the root has at least one child
                    int child = root * 2 + 1;           //root*2+1 points to the left child
                    //if the child has a sibling and the child's value is less than its sibling's...
                    if(child + 1 <= end && a[child] < a[child + 1])
                            child = child + 1;           //... then point to the right child instead
                    if(a[root] < a[child]){     //out of max-heap order
                            double tmp = a[root];
                            a[root] = a[child];
                            a[child] = tmp;
                            root = child;                //repeat to continue sifting down the child now
                    }else
                            return;
            }
    }
    
    //------------------------------------------------------------
    
    public static void heapSort(List<ProbabilisticWord> a){
	int count = a.size();
 
	//first place a in max-heap order
	heapify(a, count);
 
	int end = count - 1;
	while(end > 0){
		//swap the root(maximum value) of the heap with the
		//last element of the heap
		ProbabilisticWord tmp = a.get(end);
		a.set(end,a.get(0));
		a.set(0, tmp);
		//put the heap back in max-heap order
		siftDown(a, 0, end - 1);
		//decrement the size of the heap so that the previous
		//max value will stay in its proper place
		end--;
	}
    }

    private static void heapify(List<ProbabilisticWord> a, int count){
            //start is assigned the index in a of the last parent node
            int start = (count - 2) / 2; //binary heap

            while(start >= 0){
                    //sift down the node at index start to the proper place
                    //such that all nodes below the start index are in heap
                    //order
                    siftDown(a, start, count - 1);
                    start--;
            }
            //after sifting down the root all nodes/elements are in heap order
    }

    private static void siftDown(List<ProbabilisticWord> a, int start, int end){
            //end represents the limit of how far down the heap to sift
            int root = start;

            while((root * 2 + 1) <= end){      //While the root has at least one child
                    int child = root * 2 + 1;           //root*2+1 points to the left child
                    //if the child has a sibling and the child's value is less than its sibling's...
                    if(child + 1 <= end && a.get(child).getPpStarValue() > a.get(child+1).getPpStarValue())
                            child = child + 1;           //... then point to the right child instead
                    if(a.get(root).getPpStarValue() > a.get(child).getPpStarValue()){     //out of max-heap order
                            ProbabilisticWord tmp = a.get(root);
                            a.set(root, a.get(child));
                            a.set(child,tmp);
                            root = child;                //repeat to continue sifting down the child now
                    }else
                            return;
            }
    }
}
