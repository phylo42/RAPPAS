/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core.algos;

import java.util.*;

/**
 * This is a special implementation of SortedSet that uses a wrapper object 
 * to maintain and index that can be sorted on by a comparator, but hides 
 * that implementation detail from the external consumers so it just looks 
 * like a regular old type safe SortedSet.
 * From https://gist.github.com/jarrodhroberson/1331347
 * http://www.vertigrated.com/blog/
 * @author Jarrod Roberson
 * @param <E> 
 */
public class InsertionOrderSet<E> implements SortedSet<E>
{
    private final SortedSet<IndexedEntry<E>> bs;

    public InsertionOrderSet(final E[] e)
    {
        this(Arrays.asList(e));
    }

    public InsertionOrderSet(final List<E> l)
    {
        this();
        this.addAll(l);
    }

    public InsertionOrderSet(final Collection<? extends E> c)
    {
        this();
        this.addAll(c);
    }

    public InsertionOrderSet()
    {
        this.bs = new TreeSet<IndexedEntry<E>>(new Comparator<IndexedEntry<E>>()
        {
            @Override
            public int compare(final IndexedEntry<E> o1, final IndexedEntry<E> o2)
            {
                return o1.compareTo(o2);
            }
        });
    }

    @Override
    public boolean add(final E e)
    {
        return this.bs.add(new IndexedEntry<E>(this.bs.size(), e));
    }

    @Override
    @SuppressWarnings("unchecked")
    public boolean remove(final Object o)
    {
        return this.bs.remove(this.findByValue((E) o));
    }

    @Override
    public E first()
    {
        return this.bs.first().value;
    }

    @Override
    public E last()
    {
        return this.bs.last().value;
    }

    @Override
    public int size()
    {
        return this.bs.size();
    }

    @Override
    public Iterator<E> iterator()
    {
        return new Iterator<E>()
        {
            private Iterator<IndexedEntry<E>> it = InsertionOrderSet.this.bs.iterator();

            @Override
            public boolean hasNext()
            {
                return this.it.hasNext();
            }

            @Override
            public E next()
            {
                return this.it.next().value;
            }

            @Override
            public void remove()
            {
                this.it.remove();
            }
        };
    }

    @Override
    public boolean isEmpty()
    {
        return this.bs.isEmpty();
    }

    @Override
    public Comparator<? super E> comparator()
    {
        throw new UnsupportedOperationException("this doesn't return anything useful because all the values are wrapped with an indexer");
    }

    @Override
    @SuppressWarnings("unchecked")
    public boolean contains(final Object o)
    {
        return this.indexOf((E) o) >= 0;
    }

    @Override
    public Object[] toArray()
    {
        final Object[] a = new Object[this.bs.size()];
        for (IndexedEntry<E> ie : this.bs)
        {
            a[ie.index] = ie.value;
        }
        return a;
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> T[] toArray(final T[] a)
    {
        for (IndexedEntry<E> ie : this.bs)
        {
            a[ie.index] = (T) ie.value;
        }
        return a;
    }

    @Override
    public boolean removeAll(final Collection<?> c)
    {
        final int before = this.bs.size();
        for (Object o : c)
        {
            this.remove(o);
        }
        return this.bs.size() != before;
    }

    @Override
    public boolean addAll(final Collection<? extends E> c)
    {
        final int before = this.bs.size();
        for (E e : c)
        {
            this.add(e);
        }
        return this.bs.size() != before;
    }

    @Override
    public boolean containsAll(final Collection<?> c)
    {
        int i = 0;
        for (Object o : c)
        {
            if (this.contains(o)) { i++; }
        }
        return c.size() == i;
    }

    @Override
    public boolean retainAll(final Collection<?> c)
    {
        final int before = this.bs.size();
        final Iterator<IndexedEntry<E>> iei = this.bs.iterator();
        while (iei.hasNext())
        {
            if (!c.contains(iei.next().value)) { iei.remove(); }
        }
        return this.bs.size() != before;
    }

    @Override
    public SortedSet<E> headSet(final E toElement)
    {
        final int index = this.indexOf(toElement);
        final SortedSet<E> head = new InsertionOrderSet<E>();
        for (final IndexedEntry<E> ie : this.bs)
        {
            if (ie.index < index)
            {
                head.add(ie.value);
            }
            else
            {
                break;
            }
        }
        return head;
    }

    @Override
    public SortedSet<E> tailSet(final E fromElement)
    {
        final int index = this.indexOf(fromElement);
        final SortedSet<E> tail = new InsertionOrderSet<E>();
        for (final IndexedEntry<E> ie : this.bs)
        {
            if (ie.index >= index)
            {
                tail.add(ie.value);
            }
            else
            {
                break;
            }
        }
        return tail;
    }

    @Override
    public SortedSet<E> subSet(final E fromElement, final E toElement)
    {
        final int from = this.indexOf(fromElement);
        final int to = this.indexOf(toElement);
        final SortedSet<E> head = new InsertionOrderSet<E>();
        for (final IndexedEntry<E> ie : this.bs)
        {
            if (ie.index >= from && ie.index < to)
            {
                head.add(ie.value);
            }
            else
            {
                break;
            }
        }
        return head;
    }

    @Override
    public void clear()
    {
        this.bs.clear();
    }

    @Override
    public int hashCode()
    {
        return this.bs.hashCode();
    }

    @Override
    public boolean equals(final Object obj)
    {
        return this.bs.equals(obj);
    }

    @Override
    public String toString()
    {
        final StringBuilder sb = new StringBuilder();
        final Iterator<IndexedEntry<E>> i = this.bs.iterator();
        while (i.hasNext())
        {
            sb.append(i.next());
            if (i.hasNext()) { sb.append(","); }
        }
        return sb.toString();
    }

    private IndexedEntry<E> findByValue(final E o)
    {
        for (IndexedEntry<E> ie : this.bs)
        {
            if (ie.value.equals(o)) { return ie; }
        }
        return null;
    }

    public int indexOf(final E o)
    {
        for (IndexedEntry<E> ie : this.bs)
        {
            if (ie.value.equals(o)) { return ie.index; }
        }
        return -1;
    }

    private class IndexedEntry<E> implements Comparable<IndexedEntry<E>>
    {
        private final Integer index;
        private final E value;
        private final String strrep;

        private IndexedEntry(final int index, final E value)
        {
            this.index = index;
            this.value = value;
            this.strrep = String.format("%d:%s", this.index, this.value);
        }

        @Override
        public int compareTo(final IndexedEntry<E> o)
        {
            return this.index.compareTo(o.index);
        }

        /**
         * hashCode delegates to the contained .value hashCode implemenation
         * @return
         */
        @Override
        public int hashCode()
        {
            return this.value.hashCode();
        }

        /**
         * This returns <code>true</code> if the values are <code>equal</code>, it ignores the index
         * @param obj
         * @return
         */
        @Override
        public boolean equals(final Object obj)
        {
            return this.value.equals(obj);
        }

        @Override
        public String toString()
        {
            return this.strrep;
        }
    }
}

