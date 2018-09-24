//Seam Carving Project Code
//by Garima Kapila and Kathleen Hsu




import processing.core.PApplet;
import processing.core.PImage;

public class SeamCarving extends PApplet {
	PImage image;
	// thresholds for displaying low energy on image
	double energyThreshold = 0, energySeamThreshold = 0;

	// sets PApplet
	public void setup() {
		image = loadImage("../images/example.jpg");
		size(image.width, image.height);
		image.loadPixels();
		// thresholds initialized to be averages
		energyThreshold = avgEnergy(image);
		energySeamThreshold = avgEnergyInSeams(image);
		displaySeams(image);
	}

	// draws updated image
	public void draw() {
		image(image, 0, 0);
	}

	// displays contents of 2D double array
	public static void display(double[][] array) {
		for (int row = 0; row < array.length; row++) {
			for (int col = 0; col < array[0].length; col++) {
				System.out.print(array[row][col]);
			}
			System.out.println();
		}
	}

	// displays contents of 1D double array
	public static void display(double[] array) {
		for (int i = 0; i < array.length; i++) {
			System.out.print(array[i] + " ");
		}
	}

	// used for initializing energy threshold
	public double avgEnergy(PImage image) {
		double[][] energies = energyMatrix(image);
		double counter = 0;
		for (int r = 0; r < energies.length; r++) {
			for (int c = 0; c < energies[0].length; c++)
				counter += energies[r][c];
		}
		return counter / (energies.length * energies[0].length);
	}

	// used for initializing energy of seams threshold
	public double avgEnergyInSeams(PImage image) {
		double counter = 0;
		double[] totalEnergyInSeams = totalEnergyInSeams(image);
		for (double energy : totalEnergyInSeams) {
			counter += energy;
		}
		return counter / totalEnergyInSeams.length;
	}

	// calculates change in sub-colors between two pixels
	// uses distance formula
	public double changeInColor(int pixel1, int pixel2) {
		float red_diff = red(pixel1) - red(pixel2);
		float green_diff = green(pixel1) - green(pixel2);
		float blue_diff = blue(pixel1) - blue(pixel2);
		return red_diff * red_diff + green_diff * green_diff + blue_diff
				* blue_diff;
	}

	// calculates change in color in x- and y- gradients for pixel at (r, c)
	public double computeEnergy(int r, int c) {
		int r1 = 1, r2 = 1, c1 = 1, c2 = 1;
		if (r - 1 > 0)
			r1 = r + c * image.width - 1;
		if (r + 1 < image.width)
			r2 = r + c * image.width + 1;
		if (c - 1 > 0)
			c1 = r + (c - 1) * image.width;
		if (c + 1 < image.height)
			c2 = r + (c + 1) * image.width;
		return changeInColor(image.pixels[r1], image.pixels[r2])
				+ changeInColor(image.pixels[c1], image.pixels[c2]);
	}

	// returns array after computing energy of each pixel
	public double[][] energyMatrix(PImage image) {
		double[][] energies = new double[image.width][image.height];
		for (int r = 0; r < image.width; r++) {
			for (int c = 0; c < image.height; c++) {
				energies[r][c] = computeEnergy(r, c);
			}
		}
		return energies;
	}

	// displays spots with low energy by coloring them black
	public void lowEnergySpots(PImage image) {
		double[][] energies = energyMatrix(image);
		for (int r = 0; r < energies.length; r++) {
			for (int c = 0; c < energies[0].length; c++) {
				if (energies[r][c] < energyThreshold) {
					image.pixels[r + c * image.width] = 0;
				}
			}
		}
	}

	// returns column number of adjacent pixel with minimum energy
	public int minOfThreeAdjacent(double[][] energies, int r, int c) {
		if (r == 0)
			return c;
		if (c == 0) {
			if (Math.min(energies[r - 1][c + 1], energies[r - 1][c]) == energies[r - 1][c + 1])
				return c + 1;
			return c;
		}
		if (c == energies[0].length - 1) {
			if (Math.min(energies[r - 1][c - 1], energies[r - 1][c]) == energies[r - 1][c - 1])
				return c - 1;
			return c;
		}
		if (Math.min(Math.min(energies[r - 1][c - 1], energies[r - 1][c]),
				energies[r - 1][c + 1]) == energies[r - 1][c - 1])
			return c - 1;
		if (Math.min(Math.min(energies[r - 1][c - 1], energies[r - 1][c]),
				energies[r - 1][c + 1]) == energies[r - 1][c + 1])
			return c + 1;
		return c;
	}

	// returns array with column of adjacent pixels with min energy
	// points in the direction of where to trace seam
	public double[][] makeCostMatrix(PImage image) {
		double[][] energies = energyMatrix(image);
		double[][] costMatrix = new double[energies.length][energies[0].length];
		for (int r = 0; r < energies.length; r++) {
			for (int c = 0; c < energies[0].length; c++) {
				costMatrix[r][c] = minOfThreeAdjacent(energies, r, c);
			}
		}
		return costMatrix;
	}

	// calculates total energy within a seam
	public double totalEnergyInSeam(PImage image, int col) {
		double totalEnergy = 0;
		double[][] costMatrix = makeCostMatrix(image);
		for (int r = costMatrix.length - 1; r > 1; r--) {
			totalEnergy += computeEnergy(r - 1, (int) costMatrix[r][col]);
			col = (int) costMatrix[r][col];
		}
		return totalEnergy;
	}

	// displays seam by coloring the path in blue
	public double displaySeam(PImage image, int col) {
		double totalEnergy = totalEnergyInSeam(image, col);
		double[][] costMatrix = makeCostMatrix(image);
		if (totalEnergy < energySeamThreshold) {
			for (int r = costMatrix.length - 1; r > 1; r--) {
				image.pixels[r - 1 + ((int) (costMatrix[r][col])) * image.width] = 150;
				col = (int) costMatrix[r][col];
			}
		}
		return totalEnergy;
	}

	// calculates total energy within a seam for all possible seams
	public double[] totalEnergyInSeams(PImage image) {
		double[][] costMatrix = makeCostMatrix(image);
		double[] totalEnergyInSeams = new double[costMatrix[0].length];
		for (int c = 0; c < costMatrix[0].length - 1; c++) {
			totalEnergyInSeams[c] = totalEnergyInSeam(image, c);
		}
		return totalEnergyInSeams;
	}

	// displays all seams by coloring the paths in blue
	public double[] displaySeams(PImage image) {
		double[][] costMatrix = makeCostMatrix(image);
		double[] totalEnergyInSeams = new double[costMatrix[0].length];
		for (int c = 0; c < costMatrix[0].length - 1; c++) {
			totalEnergyInSeams[c] = totalEnergyInSeam(image, c);
			displaySeam(image, c);
		}
		return totalEnergyInSeams;
	}

	// finds pixel with lowest energy
public int findMinSeamPixel(double[][] table) {
		int lastRow = table.length - 1;
		int minIndex = 0;
		double minValue = table[lastRow][0];
		for (int j = 1; j < table[0].length; j++) {
			if (table[lastRow][j] < minValue) {
				minIndex = j;
				minValue = table[lastRow][j];
	    	}
		}
		return minValue;
	}
	// arrayList with min pixels
public Vector<Integer> findMinSeam(PImage image) {
		if (image.width == 0) {
			return new Vector<Integer>();
		}
		double[][] table = makeCostMatrix(image);
		int[] minSeamArray = new int[image.height];
		int minSeamEnd = findMinSeamPixel(table);
		minSeamArray[image.height - 1] = minSeamEnd;
		int minPath = minSeamEnd;
		for (int i = 0; i < table.length; i++) {
			minSeamArray[i] = minPath;
		minPath = findMinPath(table, i, minPath);
		}
		minSeamArray[0] = minPath;
		Vector<Integer> minSeam = new Vector<Integer>();
		for (int j = 0; j < minSeamArray.length - 1; j++) {
			minSeam.add(j);
		}
		return minSeam;
	}
	
	// finds next seam to remove
	private int findMinPath(double[][] table, int r, int c) {
		if (table[0].length == 1) {
			return 0;
		}
		if (c == 0) {
			if (table[r - 1][0] < table[r - 1][1]) {
				return 0;
			} else {
				return 1;
			}
		}
		if (c == table[0].length - 1) {
			if (table[r - 1][table[0].length - 2] < table[r - 1][table[0].length - 1])
				return table[0].length - 2;
			else
				return table[0].length - 1;
		}
		return minOfThreeAdjacent(table, r, c);
	}

//removes column of pixels
public PImage removeColumn(PImage image){
		PImage final = image.createImage(image.width - 1, image.height, PConstants.RGB);
		ArrayList<PVector> seam = findMinSeam(image);
		for(int i = 0; i < seam; i++){
			int origRowBeginning  = (int)i.y * image.width;
			int finalRowBeginning = (int)i.y * final.width; //  copy pixels on left of seam
			if(i.x > 0){
				System.arraycopy(image.pixels, origRowBeginning, result.pixels, finalRowBeginning, (int)i.x);	                                               
} //  copy pixels on right of seam
			if(i.x < image.width - 1){
				System.arraycopy(image.pixels, origRowBeginning + (int)i.x+1, final.pixels, finalRowBeginning + (int)i.x, final.width - (int)i.x);                                                      
			}
		}
		return final;
	}


}
