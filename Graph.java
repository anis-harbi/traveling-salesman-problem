package lyes;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class Graph {

  // Keep a fast index to nodes in the map
  private Map<Integer, Vertex> vertexNames;

  /**
   * Construct an empty Graph with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Graph() {
    vertexNames = new HashMap<>();
  }

  /**
   * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the graph
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the graph
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the graph
   */
  public Collection<Vertex> getVertices() {
    return vertexNames.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(String name) {
    return vertexNames.get(name);
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(int nameU, int nameV, Double cost) {
    if (!vertexNames.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertexNames.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertexNames.get(nameU);
    Vertex targetVertex = vertexNames.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param name
   *          (String) name of vertex u
   * @param name2
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(int name, int name2, double cost) {
    addEdge(name, name2, cost);
    addEdge(name2, name, cost);
  }


  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
    return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
  }

  /**
   * Computes euclidean distance between two vertices as described by their
   * coordinates
   * 
   * @param u
   *          (Vertex) vertex u
   * @param v
   *          (Vertex) vertex v
   * @return (double) distance between two vertices
   */
  public double computeEuclideanDistance(Vertex u, Vertex v) {
    return computeEuclideanDistance(u.x, u.y, v.x, v.y);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
    for (Vertex u : getVertices())
      for (Edge uv : u.adjacentEdges) {
        Vertex v = uv.target;
        uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
      }
  }



  // STUDENT CODE STARTS HERE

  
  
  
  public void generateRandomVertices(int n) {
    vertexNames = new HashMap<>(); 
    vertexNames = new HashMap<>(); 
	Random rand = new Random();
	// Create n vertices with random x, y coordinates from 0 to 100 
	for(int i = 0; i < n; i++)
	{
		int x = rand.nextInt(101);
		int y = rand.nextInt(101);
		addVertex(new Vertex(i, x, y));
	}
	// For each pair of vertices, add an undirected edge.
	for(int i = 0; i < vertexNames.size(); i++) 
		for(int j = i + 1; j < vertexNames.size(); j++) 
			addUndirectedEdge(vertexNames.get(i).name, 
					vertexNames.get(j).name, 0);
	computeAllEuclideanDistances(); // compute distances
}


  public List<Edge> nearestNeighborTsp() {
		double smallestDistance = Double.MAX_VALUE; 
		List<Edge> shortestPath = null;
		// For each vertex, perform nearest neighbor, check its length
		// against current smallest path. If smaller, replace.
		for(Vertex v : getVertices())
		{
			List<Edge> path = nearestNeighborHelper(v); 
			double length = returnPathLength(path); 
			if(length < smallestDistance)
			{
				smallestDistance = length;
				shortestPath = path;
			}
		}
		return shortestPath;
	}

	/*
	 * Helper method which performs nearest neighbor on the given vertex
	 * and returns the resulting nearest neighbor path. 
	 */
	private List<Edge> nearestNeighborHelper(Vertex v)
	{
		Vertex currentVertex = v;
		v.known = true;
		List<Edge> path = new LinkedList<>();
		
		/* get the starting vertex. Look at edges
		 * to all other unknown vertices, find shortest edge. Make the
		 * edge's target vertex known, add it to the paths list, and
		 * change currentVertex to the edge's target. Continue process
		 * with the new currentVertex.
		 */ 
		while(path.size() < vertexNames.size() - 1)
		{
			double shortestLength = Double.MAX_VALUE; 
			Edge shortestEdge = null;				 
			for(Edge e : currentVertex.adjacentEdges)
			{
				if (!e.target.known && e.distance < shortestLength)
				{
					shortestLength = e.distance;
					shortestEdge = e;
				}
			}
			shortestEdge.target.known = true;
			path.add(shortestEdge);
			currentVertex = shortestEdge.target;
		}
		for(Edge e : currentVertex.adjacentEdges)
		{
			if(e.target == v)
			{
				path.add(e);
				break;
			}
		}
		for(Vertex x : getVertices())
			x.known = false;
		return path;
	}

	/*
	 * Returns the given path's length 
	 */
	private double returnPathLength(List<Edge> path)
	{
		double sum = 0.0;
		for(Edge e : path)
			sum += e.distance;
		return sum;
	}

	
	
	
  public List<Edge> bruteForceTsp() {
		{
			// Store one permutation of vertices into array
			int[] array = new int[vertexNames.size()];
			for(int i = 0; i < array.length; i++)
				array[i] = i;
			LinkedList<int[]> permutations = new LinkedList<>();
			permute(permutations, array, 0);
			double smallestLength = Double.MAX_VALUE; // Shortest path's length
			int[] shortestPermutation = null;		  // Shortest path
			for(int[] permutation : permutations)
			{
				double totalEdgeCost = 0.0;	
				for(int i = 0; i < permutation.length; i++)
					totalEdgeCost += computeEuclideanDistance(
							vertexNames.get(permutation[i]), 
							vertexNames.get(permutation[(i+1) % 
							                permutation.length]));
				// If the total edge cost is less than current smallest, replace
				if (totalEdgeCost < smallestLength)
				{
					smallestLength = totalEdgeCost;
					shortestPermutation = permutation;
				}
			}
			
			List<Edge> shortestPath = new LinkedList<>();
			for(int i = 0; i < shortestPermutation.length; i++)
				for(Edge e:vertexNames.get(shortestPermutation[i]).adjacentEdges)
					if(e.target == vertexNames.get(shortestPermutation[(i+1) % 
					                               shortestPermutation.length]))
					{
						shortestPath.add(e);
						break;
					}
			return shortestPath;
		}

  }
	private void permute(LinkedList<int[]> permutations, int[] array, int i) 
	{
		/* We want to permute in a tree-like structure. We swap pairs of
		 * vertices and stop swapping before index i. Index i is incremented
		 * for every pair of swaps
		 */
		if (array.length - 1 == i) 
		{
			permutations.add(array.clone());
			return;
		}
				for (int j = i; j < array.length; j++) 
		{
			swap(array, i, j);
			permute(permutations, array, i+1);
			swap(array, i, j);
		}
	}
	private void swap(int[] array, int x, int y)
	{
		int temp = array[x];
		array[x] = array[y];
		array[y] = temp;
	}
  
	
	
	
	// STUDENT CODE ENDS HERE



  /**
   * Prints out the adjacency list of the graph for debugging
   */
  public void printAdjacencyList() {
    for (int u : vertexNames.keySet()) {
      StringBuilder sb = new StringBuilder();
      sb.append(u);
      sb.append(" -> [ ");
      for (Edge e : vertexNames.get(u).adjacentEdges) {
        sb.append(e.target.name);
        sb.append("(");
        sb.append(e.distance);
        sb.append(") ");
      }
      sb.append("]");
      System.out.println(sb.toString());
    }
  }
}
