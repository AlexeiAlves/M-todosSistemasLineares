package sistemas;

public class CalculoSistema {
	private CalculoSistemaAux aux;
	private double[][] ultimaI;
	private double[][] ultimaL;
	private int ultimaPerm;
	private int ultimoK;

	public CalculoSistema() {
		aux = new CalculoSistemaAux();
		ultimaPerm = 0;
		ultimaI = null;
		ultimaL = null;
		ultimoK = 0;
	}

	public double[] gauss(double[][] A, double[] b, int N) {
		double m;
		int perm = 0;

		for (int p = 0; p < N; p++) {
			perm = aux.permParcial(A, null, b, p, N);

			for (int i = p + 1; i < N; i++) {
				m = -A[i][p] / A[p][p];

				for (int j = 0; j < N; j++) {
					A[i][j] += m * A[p][j];
				}

				b[i] += m * b[p];
			}
		}

		ultimaPerm = perm;
		return aux.subsRetroativas(A, b, N);
	}

	public double[] gaussJordan(double[][] A, double[] b, int N) {
		double[][] I = aux.criarIdentidade(N);
		double m;

		for (int p = 0; p < N; p++) {
			double pivo = A[p][p];

			if (pivo != 1) {
				for (int j = 0; j < N; j++) {
					A[p][j] /= pivo;
					I[p][j] /= pivo;
				}

				b[p] /= pivo;
			}

			for (int i = 0; i < N; i++) {
				if (i != p) {
					m = -A[i][p] / A[p][p];

					for (int j = 0; j < N; j++) {
						A[i][j] += m * A[p][j];
						I[i][j] += m * I[p][j];
					}

					b[i] += m * b[p];
				}
			}
		}

		ultimaI = I;
		return b;
	}

	public double[] fatoracaoLU(double[][] A, double[] b, int N) {
		double[][] L = aux.criarIdentidade(N);
		double m;
		int perm = 0;

		for (int p = 0; p < N; p++) {
			perm = aux.permParcial(A, L, b, p, N);

			for (int i = p + 1; i < N; i++) {
				m = -A[i][p] / A[p][p];

				for (int j = 0; j < N; j++) {
					A[i][j] += m * A[p][j];
				}

				L[i][p] = -m;
			}
		}

		double[] y = aux.subsSucessivas(L, b, N);
		double[] x = aux.subsRetroativas(A, y, N);

		ultimaPerm = perm;
		ultimaL = L;
		return x;
	}

	public double[] gaussJacobi(double[][] A, double[] b, int N, double e, int maxIter) {
		double[] x0 = new double[N];
		double[] x1 = new double[N];
		double soma;

		for (int k = 0; k < maxIter; k++) {
			for (int i = 0; i < N; i++) {
				soma = 0;

				for (int j = 0; j < N && k > 0; j++) {
					if (j != i) {
						soma += A[i][j] * x0[j];
					}
				}

				x1[i] = (b[i] - soma) / A[i][i];
			}

			if (k > 0 && (aux.getNorma(x1, x0, N) < e)) {
				ultimoK = k;
				break;
			}

			x0 = x1;
			x1 = new double[N];
		}

		return x1;
	}

	public double[] gaussSeidel(double[][] A, double[] b, int N, double e, int maxIter) {
		double[] x0 = new double[N];
		double[] x1 = new double[N];
		double soma;

		for (int k = 0; k < maxIter; k++) {
			for (int i = 0; i < N; i++) {
				soma = 0;

				for (int j = 0; j < N; j++) {
					if (j != i) {
						soma += A[i][j] * x1[j];
					}
				}

				x0[i] = x1[i];
				x1[i] = (b[i] - soma) / A[i][i];
			}

			if (k > 0 && (aux.getNorma(x1, x0, N) < e)) {
				ultimoK = k;
				break;
			}
		}

		return x1;
	}

	/*---------------------------------------------------------------*/

	public double getDeterminante(double[][] U, int N) {
		return aux.getDeterminante(U, ultimaPerm, N);
	}
	
	public double[][] getUltimaInversa() {
		return ultimaI;
	}

	public double[][] getUltimaL() {
		return ultimaL;
	}

	public int getUltimoK() {
		return ultimoK;
	}
}
