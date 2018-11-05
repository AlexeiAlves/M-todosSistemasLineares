package sistemas;

public class CalculoSistemaAux {
	public double[] subsSucessivas(double[][] A, double[] b, int N) {
		double[] x = new double[N];
		double soma;

		for (int i = 0; i < N; i++) {
			soma = 0;

			for (int j = 0; j < i; j++) {
				soma += A[i][j] * x[j];
			}

			x[i] = (b[i] - soma) / A[i][i];
		}

		return x;
	}
	
	public double[] subsRetroativas(double[][] A, double[] b, int N) {
		double[] x = new double[N];
		double soma;

		for (int i = N - 1; i >= 0; i--) {
			soma = 0;

			for (int j = N - 1; j > i; j--) {
				soma += A[i][j] * x[j];
			}

			x[i] = (b[i] - soma) / A[i][i];
		}

		return x;
	}
	
	public int permParcial(double[][] U, double[][] L, double[] b, int p, int N) {
		int perm = 0, maior = p;

		for (int i = p + 1; i < N; i++) {
			maior = (Math.abs(U[i][p]) > Math.abs(U[maior][p])) ? i : maior;
		}

		if (maior != p) {
			double aux;
			perm++;

			for (int j = 0; j < N; j++) {
				aux = U[p][j];
				U[p][j] = U[maior][j];
				U[maior][j] = aux;
				
				if (L != null && j < p) {
					aux = L[p][j];
					L[p][j] = L[maior][j];
					L[maior][j] = aux;
				}
			}
			
			aux = b[p];
			b[p] = b[maior];
			b[maior] = aux;
		}
		
		return perm;
	}

	public double getNorma(double[] v, double[] x, int N) {
		double dif = 0, comp = 0, vixi, vi;

		for (int i = 0; i < N; i++) {
			vixi = Math.abs(v[i] - x[i]);
			vi = Math.abs(v[i]);

			dif = (dif < vixi) ? vixi : dif;
			comp = (comp < vi) ? vi : comp;
		}

		return dif / comp;
	}

	public double[][] criarIdentidade(int N) {
		double[][] I = new double[N][N];

		for (int i = 0; i < N; i++) {
			I[i][i] = 1;
		}

		return I;
	}

	public double[] getResiduo(double[][] A, double[] b, double[] x, int N) {
		double[] residuo = new double[N];

		for (int i = 0; i < N; i++) {
			double soma = 0;

			for (int j = 0; j < N; j++) {
				soma += (A[i][j] * x[j]);
			}

			residuo[i] = b[i] - soma;
		}

		return residuo;
	}
	
	public double getDeterminante(double[][] U, int perm, int N) {
		double det = 1;
		
		for (int p = 0; p < N; p++) {
			det *= U[p][p];
		}
		
		if ((perm % 2) == 1) {
			det *= -1;
		}
		
		return det;
	}
}
