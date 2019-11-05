package org.seqcode.projects.chexalign.alignment;

public class ArrayListSorter implements Comparable<ArrayListSorter>{
	
	public static boolean ASC = true;
	public static boolean DECS = false;
	private Object key;
	private float val;
	
	public ArrayListSorter(){}
	
	public ArrayListSorter(Object key , float val){
		super();
		this.key = key;
		this.val = val;
	}
	// Accessors
	public Object getKey(){ return key;}
	public float getVal(){ return val;}
	// Setters
	public void setKey(Object key){ this.key = key;}
	public void setVal(float val){this.val = val;}
	
	public int compareTo (ArrayListSorter compareObj){ return ((ArrayListSorter) compareObj).getVal() > this.val ? 1 : -1;} // sort in descending order
}
