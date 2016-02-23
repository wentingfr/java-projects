package prop;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class PropertyHelper {

	private static Properties properties;
	
	public static void loadProperties (String filePath) {
		//load .properties file
		InputStream in = null;
		properties = new Properties();
		
		try {
			in = new FileInputStream(new File(filePath));
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		try {
			properties.load(in);
		} 
		catch (IOException e) {
			System.out.println("Failure on property file load");
		}			
	}
	
	public static boolean HasProperty(String propertyName) {
		if(properties .containsKey(propertyName))
			return true;
		else 
			return false;
	}
	
	
	public static int getIntValue (String propertyName) throws NumberFormatException {
		return Integer.parseInt(properties.getProperty(propertyName ));			
	}
	
	public static double getDoubleValue (String propertyName) throws NumberFormatException {
		return Double.parseDouble(properties.getProperty(propertyName));			
	}
	
	public static int[] getIntArray (String propertyName, int length) throws NumberFormatException {
		int[] S = new int[length];
		String[] Sss = new String[length];
		String Ss = properties.getProperty(propertyName);

		Sss = Ss.split(";");
		
		for (int i = 0; i < length; i++) {
			S[i] = Integer.parseInt(Sss[i].trim());
		}
		
		return S;
	}
	
	public static double[] getDoubleArray (String s, int length) throws NumberFormatException {
		double[] S = new double[length];
		String[] Sss = new String[length];		
		String Ss = properties.getProperty(s);
		
		Sss = Ss.split(";");
		
		for (int i = 0; i < length; i++) {
			S[i] = Double.parseDouble(Sss[i].trim());
		}	
		
		return S;
	}
	
	public static double[][] getDouble2DArray (String propertyName, int length) throws NumberFormatException {
		double[][] e = new double[length][length];
		String[] ess = new String[length * length];		
		String es = properties.getProperty(propertyName);
		
		ess = es.split(";");
		
		for (int i = 0; i < length; i++) {
			for(int j = 0; j < length; j++) {
				e[i][j] = Double.parseDouble(ess[i * length + j].trim());
			}
		}
		
		return e;
	}
	
	public static int[][] getInt2DArray (String propertyName, int length1, int length2) throws NumberFormatException {
		int[][] S = new int[length1][length2];
		String[] Sss = new String[length1 * length2];
		String Ss = properties.getProperty(propertyName);

		Sss = Ss.split(";");
		
		for (int i = 0; i < length1; i++) {
			for(int j = 0; j < length2; j++) {
				S[i][j] = Integer.parseInt(Sss[i * length2 + j].trim());
			}
		}
		
		return S;
	}
}
