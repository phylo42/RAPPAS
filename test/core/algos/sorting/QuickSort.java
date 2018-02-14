/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos.sorting;


import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

/**
 *
 * @author ben
 */
public class QuickSort {

    /**
     * Implementation requiring Java 8
     * @param <T>
     * @param v
     * @param comparer
     * @return 
     */
    public static <T> List<T> Quicksort1(List<T> v, BiFunction<T, T, Integer> comparer) {
	if (v.size() < 2)
            return v;
	
	T pivot = v.get(v.size() / 2);

	List<T> l = new LinkedList<>(
            Quicksort1(
                v.stream().filter(
                    x -> comparer.apply(x, pivot) < 0
                ).collect(Collectors.toList()),
                comparer
            )
        );
        l.addAll(
            v.stream().filter(
                x -> comparer.apply(x, pivot) == 0
            ).collect(Collectors.toList()));
        l.addAll( 
            Quicksort1(
                v.stream().filter(
                    x -> comparer.apply(x, pivot) > 0
                ).collect(Collectors.toList()),
                comparer
            ) 
        );
	
	return l;
    }
    
    /**
     * more general implementation
     * @param numbers
     * @return 
     */
    public static ArrayList<Double> quicksort2(ArrayList<Double> numbers) {
        if (numbers.size() <= 1)
            return numbers;
        int pivot = numbers.size() / 2;
        ArrayList<Double> lesser = new ArrayList<Double>();
        ArrayList<Double> greater = new ArrayList<Double>();
        int sameAsPivot = 0;
        for (double number : numbers) {
            if (number > numbers.get(pivot))
                greater.add(number);
            else if (number < numbers.get(pivot))
                lesser.add(number);
            else
                sameAsPivot++;
        }
        lesser = quicksort2(lesser);
        for (int i = 0; i < sameAsPivot; i++)
            lesser.add(numbers.get(pivot));
        greater = quicksort2(greater);
        ArrayList<Double> sorted = new ArrayList<Double>();
        for (double number : lesser)
            sorted.add(number);
        for (double number: greater)
            sorted.add(number);
        return sorted;
    }
    
    //#######################################
    //below implementation using a randomly selected pivot
    public static final Random RND = new Random();

    private static void swap(Object[] array, int i, int j) {
        Object tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }

    private static <E extends Comparable<? super E>> int partition(E[] array, int begin, int end) {
        int index = begin + RND.nextInt(end - begin + 1);
        E pivot = array[index];
        swap(array, index, end);
        for (int i = index = begin; i < end; ++i) {
            if (array[i].compareTo(pivot) <= 0) {
                swap(array, index++, i);
            }
        }
        swap(array, index, end);
        return (index);
    }

    private static <E extends Comparable<? super E>> void qsort(E[] array, int begin, int end) {
        if (end > begin) {
            int index = partition(array, begin, end);
            qsort(array, begin, index - 1);
            qsort(array, index + 1, end);
        }
    }

    public static <E extends Comparable<? super E>> void quickSort3(E[] array) {
        qsort(array, 0, array.length - 1);
    }
    
    @SuppressWarnings("unchecked")
    private static <E extends Comparable<? super E>> List<E>[] split(List<E> v) {
            List<E>[] results = (List<E>[])new List[] { new LinkedList<E>(), new LinkedList<E>() };
            Iterator<E> it = v.iterator();
            E pivot = it.next();
            while (it.hasNext()) {
                    E x = it.next();
                    if (x.compareTo(pivot) < 0) results[0].add(x);
                    else results[1].add(x);
            }
            return results;
    }

    public static <E extends Comparable<? super E>> List<E> quicksort4(List<E> v) {
            if (v == null || v.size() <= 1) return v;
            final List<E> result = new LinkedList<E>();
            final List<E>[] lists = split(v);
            result.addAll(quicksort4(lists[0]));
            result.add(v.get(0));
            result.addAll(quicksort4(lists[1]));
            return result;
    }
    
}
