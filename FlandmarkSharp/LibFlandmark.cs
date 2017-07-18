using OpenCvSharp;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlandmarkSharp
{
    class FLANDMARK_PSIG
    {
        public int[] disp;
        public int ROWS, COLS;
    }

    class FLANDMARK_Options
    {
        public byte M;
        public int[] S;
        public int[] bw = new int[2];
        public int[] bw_margin = new int[2];
        public FLANDMARK_PSIG[] PsiGS0, PsiGS1, PsiGS2;
        public int[] PSIG_ROWS = new int[3];
        public int[] PSIG_COLS = new int[3];
    }

    class FLANDMARK_LBP
    {
        public int[] winSize = new int[2];
        public byte hop;
        public uint[] wins;
        public uint WINS_ROWS, WINS_COLS;
    }

    class FLANDMARK_Data
    {
        public FLANDMARK_LBP[] lbp;
        public int[] imSize = new int[2];
        public int[] mapTable;
        public FLANDMARK_Options options = new FLANDMARK_Options();
    }

    class FLANDMARK_Model
    {
        public double[] W;
        public int W_ROWS, W_COLS;
        public FLANDMARK_Data data = new FLANDMARK_Data();
        public byte[] normalizedImageFrame;
        public double[] bb;
        public float[] sf;
    }

    class FLANDMARK_PSI
    {
        public char[] data;
        public uint PSI_ROWS, PSI_COLS;
    }

    class FLANDMARK_PSI_SPARSE
    {
        public uint[] idxs;
        public uint PSI_ROWS, PSI_COLS;
    }

    public enum EError_T
    {
        NO_ERR = 0,
        ERROR_M = 1,
        ERROR_BW = 2,
        ERROR_BW_MARGIN = 3,
        ERROR_W = 4,
        ERROR_DATA_IMAGES = 5,
        ERROR_DATA_MAPTABLE = 6,
        ERROR_DATA_LBP = 7,
        ERROR_DATA_OPTIONS_S = 8,
        ERROR_DATA_OPTIONS_PSIG = 9,
        UNKNOWN_ERROR = 100
    };

    public class LibFlandmark
    {
        FLANDMARK_Model model;

        static int DefINDEX(int r, int c, int nr)
        {
            return c * nr + r;
        }

        static int DefROW(int i, int r)
        {
            return (i - 1) % r;
        }

        static int DefCOL(int i, int r)
        {
            return (i - 1) / r;
        }

        private string readUntilSpace(StreamReader reader)
        {
            char current = 'A';
            StringBuilder builder = new StringBuilder();

            while (current != ' ')
            {
                current = (char)reader.Read();
                builder.Append(current);
            }

            string ret = builder.ToString();
            Console.WriteLine(ret);
            return ret;
        }

        public LibFlandmark(string filename)
        {
            int[] p_int;
            int tsize = -1, tmp_tsize = -1;
            int[] p_uint8;

            // TODO: using
            FileStream fin = File.Open(filename, FileMode.Open, FileAccess.ReadWrite);
            StreamReader reader = new StreamReader(fin, Encoding.ASCII);

            // allocate memory for FLANDMARK_Model
            FLANDMARK_Model tst = new FLANDMARK_Model();

            fin.Position = 0;

            readUntilSpace(reader);
            tst.data.options.M = (byte)readUntilSpace(reader)[0];

            readUntilSpace(reader);
            tst.data.options.bw[0] = Convert.ToInt32(readUntilSpace(reader));
            tst.data.options.bw[1] = Convert.ToInt32(readUntilSpace(reader));

            readUntilSpace(reader);
            tst.data.options.bw_margin[0] = Convert.ToInt32(readUntilSpace(reader));
            tst.data.options.bw_margin[1] = Convert.ToInt32(readUntilSpace(reader));

            readUntilSpace(reader);
            tst.W_ROWS = Convert.ToInt32(readUntilSpace(reader));
            tst.W_COLS = Convert.ToInt32(readUntilSpace(reader));

            readUntilSpace(reader);
            tst.data.imSize[0] = Convert.ToInt32(readUntilSpace(reader));
            tst.data.imSize[1] = Convert.ToInt32(readUntilSpace(reader));

            int M = tst.data.options.M;

            tst.data.lbp = new FLANDMARK_LBP[M];

            for (int i = 0; i < M; i++)
            {
                readUntilSpace(reader);

                var lbp = new FLANDMARK_LBP();

                lbp.WINS_ROWS = (uint)Convert.ToInt32(readUntilSpace(reader));
                lbp.WINS_COLS = (uint)Convert.ToInt32(readUntilSpace(reader));

                tst.data.lbp[i] = lbp;
            }

            for (int i = 0; i < 3; i++)
            {
                readUntilSpace(reader);
                tst.data.options.PSIG_ROWS[i] = Convert.ToInt32(readUntilSpace(reader));
                tst.data.options.PSIG_COLS[i] = Convert.ToInt32(readUntilSpace(reader));
            }

            fin.Close();
            fin.Dispose();

            fin = File.Open(filename, FileMode.Open, FileAccess.ReadWrite);
            fin.Position = 111;
            BinaryReader breader = new BinaryReader(fin);

            //            // load model.W -----------------------------------------------------------
            tst.W = new double[tst.W_ROWS];

            for (int i = 0; i < tst.W_ROWS; i++)
            {
                tst.W[i] = breader.ReadDouble();
            }

            //            // load model.data.mapTable -----------------------------------------------
            p_int = new int[M * 4];
            tst.data.mapTable = new int[M * 4];

            for (int i = 0; i < M * 4; i++)
            {
                p_int[i] = breader.ReadInt32();
            }

            //            for (int i = 0; i < M * 4; ++i)
            //            {
            //                tst->data.mapTable[i] = p_int[i];
            //            }
            for (int i = 0; i < M * 4; i++)
            {
                tst.data.mapTable[i] = p_int[i];
            }
            p_int = null;

            //            // load model.data.lbp ---------------------------------------------------
            for (int i = 0; i < M; i++)
            {
                //                // lbp{idx}.winSize
                //                p_int = (int*)malloc(2 * sizeof(int));
                p_int = new int[2];

                //                if (fread(p_int, 2 * sizeof(int), 1, fin) != 1)
                //                {
                //                    Console.Write("Error reading file %s\n", filename);
                //                    return 0;
                //                }
                p_int[0] = breader.ReadInt32();
                p_int[1] = breader.ReadInt32();

                //                for (int i = 0; i < 2; ++i)
                //                {
                //                    tst->data.lbp[idx].winSize[i] = p_int[i];
                //                }
                tst.data.lbp[i].winSize[0] = p_int[0];
                tst.data.lbp[i].winSize[1] = p_int[1];

                //                free(p_int);
                p_int = null;

                //                // lbp{idx}.hop
                //                p_uint8 = (uint8_t*)malloc(sizeof(uint8_t));
                //                if (fread(p_uint8, sizeof(uint8_t), 1, fin) != 1)
                //                {
                //                    Console.Write("Error reading file %s\n", filename);
                //                    return 0;
                //                }
                //                tst->data.lbp[idx].hop = p_uint8[0];
                //                free(p_uint8);
                tst.data.lbp[i].hop = breader.ReadByte();

                //                // lbp{idx}.wins
                //                tsize = tst->data.lbp[idx].WINS_ROWS * tst->data.lbp[idx].WINS_COLS;
                tsize = (int)(tst.data.lbp[i].WINS_ROWS * tst.data.lbp[i].WINS_COLS);

                //                tst->data.lbp[idx].wins = (uint32_t*)malloc(tsize * sizeof(uint32_t));
                tst.data.lbp[i].wins = new uint[tsize];

                //                if (fread(tst->data.lbp[idx].wins, tsize * sizeof(uint32_t), 1, fin) != 1)
                //                {
                //                    Console.Write("Error reading file %s\n", filename);
                //                    return 0;
                //                    //exit(1);
                //                }
                for (int r = 0; r < tsize; r++)
                {
                    tst.data.lbp[i].wins[r] = breader.ReadUInt32();
                }
            }

            //            // load model.options.S --------------------------------------------------
            tst.data.options.S = new int[4 * M];

            //            if (fread(p_int, 4 * M * sizeof(int), 1, fin) != 1)
            //            {
            //                Console.Write("Error reading file %s\n", filename);
            //                return 0;
            //            }
            //            for (int i = 0; i < 4 * M; ++i)
            //            {
            //                tst->data.options.S[i] = p_int[i];
            //            }
            for (int i = 0; i < 4 * M; i++)
            {
                tst.data.options.S[i] = breader.ReadInt32();
            }
            p_int = null;

            //            // load model.options.PsiG -----------------------------------------------
            FLANDMARK_PSIG[] PsiGi = null;

            for (int psigs_ind = 0; psigs_ind < 3; psigs_ind++)
            {
                tsize = tst.data.options.PSIG_ROWS[psigs_ind] * tst.data.options.PSIG_COLS[psigs_ind];

                switch (psigs_ind)
                {
                    case 0:
                        //                        tst->data.options.PsiGS0 = (FLANDMARK_PSIG*)malloc(tsize * sizeof(FLANDMARK_PSIG));
                        //                        PsiGi = tst->data.options.PsiGS0;
                        tst.data.options.PsiGS0 = new FLANDMARK_PSIG[tsize];
                        PsiGi = tst.data.options.PsiGS0;
                        break;
                    case 1:
                        //                        tst->data.options.PsiGS1 = (FLANDMARK_PSIG*)malloc(tsize * sizeof(FLANDMARK_PSIG));
                        //                        PsiGi = tst->data.options.PsiGS1;
                        tst.data.options.PsiGS1 = new FLANDMARK_PSIG[tsize];
                        PsiGi = tst.data.options.PsiGS1;
                        break;
                    case 2:
                        //                        tst->data.options.PsiGS2 = (FLANDMARK_PSIG*)malloc(tsize * sizeof(FLANDMARK_PSIG));
                        //                        PsiGi = tst->data.options.PsiGS2;
                        tst.data.options.PsiGS2 = new FLANDMARK_PSIG[tsize];
                        PsiGi = tst.data.options.PsiGS2;
                        break;
                        //                }
                }

                //                for (int idx = 0; idx < tsize; ++idx)
                //                {
                int temp = 0;
                for (int i = 0; i < tsize; i++)
                {
                    PsiGi[i] = new FLANDMARK_PSIG();

                    //                    // disp ROWS
                    //                    p_int = (int*)malloc(sizeof(int));
                    //                    if (fread(p_int, sizeof(int), 1, fin) != 1)
                    //                    {
                    //                        Console.Write("Error reading file %s\n", filename);
                    //                        return 0;
                    //                        //exit(1);
                    //                    }
                    temp = breader.ReadInt32();

                    //                    PsiGi[idx].ROWS = p_int[0];
                    PsiGi[i].ROWS = temp;
                    //                    free(p_int);

                    //                    // disp COLS
                    //                    p_int = (int*)malloc(sizeof(int));
                    //                    if (fread(p_int, sizeof(int), 1, fin) != 1)
                    //                    {
                    //                        Console.Write("Error reading file %s\n", filename);
                    //                        return 0;
                    //                        //exit(1);
                    //                    }
                    temp = breader.ReadInt32();

                    //                    PsiGi[idx].COLS = p_int[0];
                    PsiGi[i].COLS = temp;
                    //                    free(p_int);

                    //                    // disp
                    //                    tmp_tsize = PsiGi[idx].ROWS * PsiGi[idx].COLS;
                    tmp_tsize = PsiGi[i].ROWS * PsiGi[i].COLS;

                    //                    PsiGi[idx].disp = (int*)malloc(tmp_tsize * sizeof(int));
                    PsiGi[i].disp = new int[tmp_tsize];

                    //                    if (fread(PsiGi[idx].disp, tmp_tsize * sizeof(int), 1, fin) != 1)
                    //                    {
                    //                        Console.WriteLine("Error reading file %s\n", filename);
                    //                        return 0;
                    //                        //exit(1);
                    //                    }
                    for (int r = 0; r < tmp_tsize; r++)
                    {
                        PsiGi[i].disp[r] = breader.ReadInt32();
                    }
                }
            }

            //            // TODO: using
            fin.Close();

            //            tst->normalizedImageFrame = (uint8_t*)calloc(tst->data.options.bw[0] * tst->data.options.bw[1], sizeof(uint8_t));
            tst.normalizedImageFrame = new byte[tst.data.options.bw[0] * tst.data.options.bw[1]];

            //            tst.bb = (double*)calloc(4, sizeof(double));
            tst.bb = new double[4];

            //            tst.sf = (float*)calloc(2, sizeof(float));
            tst.sf = new float[2];

            model = tst;
        }

        // NEED
        //void flandmark_get_psi_mat_sparse(FLANDMARK_PSI_SPARSE* Psi, FLANDMARK_Model* model, int lbpidx)
        //{
        void get_psi_mat_sparse(ref FLANDMARK_PSI_SPARSE Psi, ref FLANDMARK_Model model, int lbpidx)
        {
            //    t_index* Features;
            //    uint8_t* Images = model->normalizedImageFrame;
            //    uint32_t im_H = (uint32_t)model->data.imSize[0];
            //    uint32_t im_W = (uint32_t)model->data.imSize[1];
            //    uint32_t* Wins = model->data.lbp[lbpidx].wins;
            //    uint16_t win_H = (uint16_t)model->data.lbp[lbpidx].winSize[0];
            //    uint16_t win_W = (uint16_t)model->data.lbp[lbpidx].winSize[1];
            //    uint16_t nPyramids = model->data.lbp[lbpidx].hop;
            //    uint32_t nDim = liblbp_pyr_get_dim(win_H, win_W, nPyramids) / 256;
            //    uint32_t nData = model->data.lbp[lbpidx].WINS_COLS;
            uint[] Features;
            byte[] Images = model.normalizedImageFrame;
            uint im_H = (uint)model.data.imSize[0];
            uint im_W = (uint)model.data.imSize[1];
            uint[] Wins = model.data.lbp[lbpidx].wins;
            UInt16 win_H = (UInt16)model.data.lbp[lbpidx].winSize[0];
            UInt16 win_W = (UInt16)model.data.lbp[lbpidx].winSize[1];
            UInt16 nPyramids = model.data.lbp[lbpidx].hop;
            uint nDim = LibLBP.PyrGetDim(win_H, win_W, nPyramids) / 256;
            uint nData = model.data.lbp[lbpidx].WINS_COLS;

            //    uint32_t cnt0, mirror, x, x1, y, y1, idx, *win;
            //    uint8_t* img_ptr;
            uint cnt0, mirror, x, x1, y, y1, idx;
            uint[] win;
            // TODO
            uint img_ptr;

            //    Features = (t_index*)calloc(nDim * nData, sizeof(t_index));
            Features = new uint[nDim * nData];
            //    if (Features == NULL)
            //    {
            //        Console.Write("Not enough memory for LBP features.\n");
            //        exit(1);
            //    }

            //    win = (uint32_t*)calloc(win_H * win_W, sizeof(uint32_t));
            win = new uint[win_H * win_W];
            //    if (win == NULL)
            //    {
            //        Console.Write("Not enough memory for cropped_window.\n");
            //        exit(1);
            //    }

            //    for (uint32_t i = 0; i < nData; ++i)
            //    {
            for (uint i = 0; i < nData; i++)
            {
                //        idx = Wins[INDEX(0, i, 4)] - 1;
                //        x1 = Wins[INDEX(1, i, 4)] - 1;
                //        y1 = Wins[INDEX(2, i, 4)] - 1;
                //        mirror = Wins[INDEX(3, i, 4)];
                idx = Wins[DefINDEX(0, (int)i, 4)] - 1;
                x1 = Wins[DefINDEX(1, (int)i, 4)] - 1;
                y1 = Wins[DefINDEX(2, (int)i, 4)] - 1;
                mirror = Wins[DefINDEX(3, (int)i, 4)];

                //        img_ptr = &Images[idx * im_H * im_W];
                img_ptr = idx * im_H * im_W;

                //        cnt0 = 0;
                cnt0 = 0;

                //        if (mirror == 0)
                //        {
                //            for (x = x1; x < x1 + win_W; x++)
                //                for (y = y1; y < y1 + win_H; y++)
                //                    win[cnt0++] = img_ptr[INDEX(y, x, im_H)];
                //        }
                //        else
                //        {
                //            for (x = x1 + win_W - 1; x >= x1; x--)
                //                for (y = y1; y < y1 + win_H; y++)
                //                    win[cnt0++] = img_ptr[INDEX(y, x, im_H)];
                //        }
                if (mirror == 0)
                {
                    for (x = x1; x < x1 + win_W; x++)
                        for (y = y1; y < y1 + win_H; y++)
                            win[cnt0++] = Images[img_ptr + DefINDEX((int)y, (int)x, (int)im_H)];
                }
                else
                {
                    for (x = x1 + win_W - 1; x >= x1; x--)
                        for (y = y1; y < y1 + win_H; y++)
                            win[cnt0++] = Images[img_ptr + DefINDEX((int)y, (int)x, (int)im_H)];
                }

                //        liblbp_pyr_features_sparse(&Features[nDim * i], nDim, win, win_H, win_W);
                LibLBP.PyrFeaturesSparse(Features, nDim, win, win_H, win_W, nDim * i);
                //    }
            }

            //    Psi->PSI_COLS = nData;
            //    Psi->PSI_ROWS = nDim;
            //    Psi->idxs = Features;
            Psi.PSI_COLS = nData;
            Psi.PSI_ROWS = nDim;
            Psi.idxs = Features;

            //    free(win);
            //}
        }

        // NEED
        //int flandmark_detect_base(uint8_t* face_image, FLANDMARK_Model* model, double* landmarks)
        //{
        void detect_base(ref byte[] face_image, ref FLANDMARK_Model model, ref double[] landmarks)
        {
            //    const int M = model->data.options.M;
            int M = model.data.options.M;

            //    const double* W = model->W;
            double[] W = model.W;

            //    int tsize = -1, cols = -1, rows = -1;
            int tsize = -1, cols = -1, rows = -1;

            //    const int* mapTable = model->data.mapTable;
            int[] mapTable = model.data.mapTable;

            //    if (!model->normalizedImageFrame)
            //    {
            //        model->normalizedImageFrame = face_image;
            //    }
            //    else
            //    {
            //        //
            //    }
            if (model.normalizedImageFrame == null)
            {
                model.normalizedImageFrame = face_image;
            }
            else
            {
                //
            }

            //    // get PSI matrix
            //    FLANDMARK_PSI_SPARSE* Psi_sparse = (FLANDMARK_PSI_SPARSE*)malloc(M * sizeof(FLANDMARK_PSI_SPARSE));
            FLANDMARK_PSI_SPARSE[] Psi_sparse = new FLANDMARK_PSI_SPARSE[M];

            //    for (int idx = 0; idx < M; ++idx)
            //    {
            //        flandmark_get_psi_mat_sparse(&Psi_sparse[idx], model, idx);
            //    }
            for (int i = 0; i < M; i++)
            {
                Psi_sparse[i] = new FLANDMARK_PSI_SPARSE();
                get_psi_mat_sparse(ref Psi_sparse[i], ref model, i);
            }

            //    // get Q and G
            //    double** q = (double**)calloc(M, sizeof(double*));
            //    double** g = (double**)calloc((M - 1), sizeof(double*));
            List<double[]> q = new List<double[]>(M);
            List<double[]> g = new List<double[]>(M - 1);

            //    int idx_qtemp = 0;
            int idx_qtemp = 0;

            //    for (int idx = 0; idx < M; ++idx)
            //    {
            for (int i = 0; i < M; i++)
            {
                //        // Q
                //        tsize = mapTable[INDEX(idx, 1, M)] - mapTable[INDEX(idx, 0, M)] + 1;
                tsize = mapTable[DefINDEX(i, 1, M)] - mapTable[DefINDEX(i, 0, M)] + 1;

                //        double* q_temp = (double*)calloc(tsize, sizeof(double));
                double[] q_temp = new double[tsize];
                //        memcpy(q_temp, W + mapTable[INDEX(idx, 0, M)] - 1, tsize * sizeof(double));
                Array.Copy(W, mapTable[DefINDEX(i, 0, M)], q_temp, 0, tsize);

                //        // sparse dot product <W_q, PSI_q>
                //        cols = Psi_sparse[idx].PSI_COLS; rows = Psi_sparse[idx].PSI_ROWS;
                cols = (int)Psi_sparse[i].PSI_COLS;
                rows = (int)Psi_sparse[i].PSI_ROWS;
                //        uint32_t* psi_temp = Psi_sparse[idx].idxs;
                uint[] psi_temp = Psi_sparse[i].idxs;

                //        q[idx] = (double*)malloc(cols * sizeof(double));
                q.Add(new double[cols]);

                //        for (int i = 0; i < cols; ++i)
                //        {
                //            double dotprod = 0.0f;
                //            for (int j = 0; j < rows; ++j)
                //            {
                //                idx_qtemp = psi_temp[(rows * i) + j];
                //                dotprod += q_temp[idx_qtemp];
                //            }
                //            q[idx][i] = dotprod;
                //        }
                //        free(q_temp);
                for (int c = 0; c < cols; c++)
                {
                    double dotprod = 0.0f;
                    for (int j = 0; j < rows; ++j)
                    {
                        idx_qtemp = (int)psi_temp[(rows * c) + j];
                        dotprod += q_temp[idx_qtemp];
                    }
                    q[i][c] = dotprod;
                }

                //        // G
                //        if (idx > 0)
                //        {
                //            tsize = mapTable[INDEX(idx, 3, M)] - mapTable[INDEX(idx, 2, M)] + 1;
                //            g[idx - 1] = (double*)malloc(tsize * sizeof(double));
                //            memcpy(g[idx - 1], W + mapTable[INDEX(idx, 2, M)] - 1, tsize * sizeof(double));
                //        }
                if (i > 0)
                {
                    tsize = mapTable[DefINDEX(i, 3, M)] - mapTable[DefINDEX(i, 2, M)] + 1;
                    g.Add(new double[1]);
                    g[i - 1] = new double[tsize];
                    Array.Copy(W, mapTable[DefINDEX(i, 2, M)] - 1, g[i - 1], 0, tsize);
                }

                //    }
            }

            //    // argmax
            //    flandmark_argmax(landmarks, &model->data.options, mapTable, Psi_sparse, q, g);
            argmax(landmarks, ref model.data.options, mapTable, Psi_sparse, q, g);

            //    // cleanup Psi_sparse[].idxs
            //    for (int i = 0; i < M; ++i)
            //    {
            //        free(Psi_sparse[i].idxs);
            //    }
            //    free(Psi_sparse);

            //    // cleanup q
            //    for (int i = 0; i < M; ++i)
            //    {
            //        free(q[i]);
            //    }
            //    free(q);
            //    // cleanup g
            //    for (int i = 0; i < M - 1; ++i)
            //    {
            //        free(g[i]);
            //    }
            //    free(g);

            //    return 0;
            //}
        }

        // NEED
        //int flandmark_detect(IplImage* img, int* bbox, FLANDMARK_Model* model, double* landmarks, int* bw_margin)
        //{
        public double[] Detect(Mat img, int[] bbox, int?[] bw_margin = null)
        {
            double[] landmarks = new double[model.data.options.M * 2];
            //    int retval = 0;

            //    if (bw_margin)
            //    {
            //        model->data.options.bw_margin[0] = bw_margin[0];
            //        model->data.options.bw_margin[1] = bw_margin[1];
            //    }
            if (bw_margin != null)
            {
                model.data.options.bw_margin[0] = (int)bw_margin[0];
                model.data.options.bw_margin[1] = (int)bw_margin[1];
            }

            //    // Get normalized image frame
            //    retval = flandmark_get_normalized_image_frame(img, bbox, model->bb, model->normalizedImageFrame, model);
            get_normalized_image_frame(img, bbox, model.bb, model.normalizedImageFrame, model);
            //    if (retval)
            //    {
            //        // flandmark_get_normlalized_image_frame ERROR;
            //        return 1;
            //    }

            //    // Call flandmark_detect_base
            //    retval = flandmark_detect_base(model->normalizedImageFrame, model, landmarks);
            detect_base(ref model.normalizedImageFrame, ref model, ref landmarks);
            //    if (retval)
            //    {
            //        // flandmark_detect_base ERROR
            //        return 2;
            //    }

            //    // transform coordinates of detected landmarks from normalized image frame back to the original image
            //    model->sf[0] = (float)(model->bb[2] - model->bb[0]) / model->data.options.bw[0];
            //    model->sf[1] = (float)(model->bb[3] - model->bb[1]) / model->data.options.bw[1];
            model.sf[0] = (float)(model.bb[2] - model.bb[0]) / model.data.options.bw[0];
            model.sf[1] = (float)(model.bb[3] - model.bb[1]) / model.data.options.bw[1];
            //    for (int i = 0; i < 2 * model->data.options.M; i += 2)
            //    {
            //        landmarks[i] = landmarks[i] * model->sf[0] + model->bb[0];
            //        landmarks[i + 1] = landmarks[i + 1] * model->sf[1] + model->bb[1];
            //    }
            for (int i = 0; i < 2 * model.data.options.M; i += 2)
            {
                landmarks[i] = landmarks[i] * model.sf[0] + model.bb[0];
                landmarks[i + 1] = landmarks[i + 1] * model.sf[1] + model.bb[1];
            }

            //    return 0;
            return landmarks;

            //}
        }

        //// NEED
        //////int flandmark_imcrop(IplImage *input, IplImage *output, const CvRect region, FLANDMARK_Model *model)
        ////int flandmark_imcrop(IplImage* input, IplImage* output, const CvRect region)
        ////{
        //void imcrop(Mat input, Mat output, Rect region)
        //{
        //    //    if (input->width <= 0 || input->height <= 0 || region.width <= 0 || region.height <= 0)
        //    //    {
        //    //        return 1;
        //    //    }
        //    if (input.Width <= 0 || input.Height <= 0 || region.Width <= 0 || region.Height <= 0)
        //    {
        //        throw new Exception();
        //    }

        //    //    if (input->depth != IPL_DEPTH_8U)
        //    //    {
        //    //        return 1;
        //    //    }

        //    //    cvSetImageROI(input, region);

        //    //    if (output->width < region.width || output->height < region.height)
        //    //    {
        //    //        cvReleaseImage(&output);
        //    //        output = cvCreateImage(cvSize(region.width, region.height), IPL_DEPTH_8U, 1);
        //    //    }
        //    //    else
        //    //    {
        //    //        output->width = region.width;
        //    //        output->height = region.height;
        //    //    }
        //    if(output.Width < region.Width || output.Height < region.Height)
        //    {
        //        output = new Mat(region.Size, MatType.CV_8UC1);
        //    }
        //    else
        //    {
        //        output.Resize(region.Size);
        //    }

        //    //    cvCopy(input, output, NULL);
        //    //    cvResetImageROI(input);

        //    //    return 0;
        //    //}
        //}

        // NEED
        //int flandmark_get_normalized_image_frame(IplImage* input, const int bbox[], double* bb, uint8_t* face_img, FLANDMARK_Model* model)
        //{
        void get_normalized_image_frame(Mat input, int[] bbox, double[] bb, byte[] face_img, FLANDMARK_Model model)
        {
            //    bool flag;
            bool flag;
            //    int d[2];
            int[] d = new int[2];
            //    double c[2], nd[2];
            double[] c = new double[2], nd = new double[2];

            //    // extend bbox by bw_margin
            //    d[0] = bbox[2] - bbox[0] + 1; d[1] = bbox[3] - bbox[1] + 1;
            //    c[0] = (bbox[2] + bbox[0]) / 2.0f; c[1] = (bbox[3] + bbox[1]) / 2.0f;
            //    nd[0] = d[0] * model->data.options.bw_margin[0] / 100.0f + d[0];
            //    nd[1] = d[1] * model->data.options.bw_margin[1] / 100.0f + d[1];
            d[0] = bbox[2] - bbox[0] + 1;
            d[1] = bbox[3] - bbox[1] + 1;
            c[0] = (bbox[2] + bbox[0]) / 2.0f;
            c[1] = (bbox[3] + bbox[1]) / 2.0f;
            nd[0] = d[0] * model.data.options.bw_margin[0] / 100.0f + d[0];
            nd[1] = d[1] * model.data.options.bw_margin[1] / 100.0f + d[1];

            //    bb[0] = (c[0] - nd[0] / 2.0f);
            //    bb[1] = (c[1] - nd[1] / 2.0f);
            //    bb[2] = (c[0] + nd[0] / 2.0f);
            //    bb[3] = (c[1] + nd[1] / 2.0f);
            bb[0] = (c[0] - nd[0] / 2.0f);
            bb[1] = (c[1] - nd[1] / 2.0f);
            bb[2] = (c[0] + nd[0] / 2.0f);
            bb[3] = (c[1] + nd[1] / 2.0f);

            //    flag = bb[0] > 0 && bb[1] > 0 && bb[2] < input->width && bb[3] < input->height
            //        && bbox[0] > 0 && bbox[1] > 0 && bbox[2] < input->width && bbox[3] < input->height;
            flag = bb[0] > 0 && bb[1] > 0 && bb[2] < input.Width && bb[3] < input.Height && bbox[0] > 0 && bbox[1] > 0 && bbox[2] < input.Width && bbox[3] < input.Height;

            //    if (!flag)
            //    {
            //        return 1;
            //    }
            if (!flag)
                throw new Exception();

            //    IplImage* croppedImage = cvCreateImage(cvSize(input->width, input->height), IPL_DEPTH_8U, 1);
            //    IplImage* resizedImage = cvCreateImage(cvSize(model->data.options.bw[0], model->data.options.bw[1]), IPL_DEPTH_8U, 1);

            //    // crop and resize image to normalized frame
            //    if (flandmark_imcrop(input, croppedImage, cvRect((int)bb[0], (int)bb[1], (int)bb[2] - (int)bb[0] + 1, (int)bb[3] - (int)bb[1] + 1)))
            //    {
            //        // something was bad
            //        return 1;
            //  
            Mat resizedImage = new Mat(input, new Rect((int)bb[0], (int)bb[1], (int)bb[2] - (int)bb[0] + 1, (int)bb[3] - (int)bb[1] + 1));
            resizedImage = resizedImage.Resize(new Size(model.data.options.bw[0], model.data.options.bw[1]));

            //    // resize
            //    cvResize(croppedImage, resizedImage, CV_INTER_CUBIC);

            //    // tranform IplImage to simple 1D uint8 array representing 2D uint8 normalized image frame
            //    for (int x = 0; x < model->data.options.bw[0]; ++x)
            //    {
            for (int x = 0; x < model.data.options.bw[0]; x++)
            {
                //        for (int y = 0; y < model->data.options.bw[1]; ++y)
                //        {
                //            face_img[INDEX(y, x, model->data.options.bw[1])] = (uint8_t)((resizedImage->imageData + resizedImage->widthStep * y)[x]);
                //        }
                for (int y = 0; y < model.data.options.bw[1]; y++)
                {
                    face_img[DefINDEX(x, y, model.data.options.bw[1])] = resizedImage.At<byte>(x, y);
                }

                //    }
            }

            //    cvReleaseImage(&croppedImage);
            //    cvReleaseImage(&resizedImage);
            // TODO: using
            resizedImage.Dispose();

            //    return 0;
            //}
        }


        //void flandmark_maximize_gdotprod(double* maximum, double* idx, const double* first, const double* second, const int* third, const int cols, const int tsize)
        //{
        //    *maximum = -FLT_MAX;
        //    *idx = -1;
        //    for (int dp_i = 0; dp_i < cols; ++dp_i)
        //    {
        //        double dotprod = 0.0f;
        //        for (int dp_j = 0; dp_j < tsize; ++dp_j)
        //        {
        //            dotprod += second[dp_j] * (double)(third[dp_i * tsize + dp_j]);
        //        }
        //        if (*maximum < first[dp_i] + dotprod)
        //        {
        //            *idx = dp_i;
        //            *maximum = first[dp_i] + dotprod;
        //        }
        //    }
        //}
        void maximize_gdotprod(ref double maximum, ref double idx, double[] first, double[] second, int[] third, int cols, int tsize)
        {
            maximum = -float.MaxValue;
            idx = -1;
            for (int dp_i = 0; dp_i < cols; ++dp_i)
            {
                double dotprod = 0.0f;
                for (int dp_j = 0; dp_j < tsize; ++dp_j)
                {
                    dotprod += second[dp_j] * (double)(third[dp_i * tsize + dp_j]);
                }
                if (maximum < first[dp_i] + dotprod)
                {
                    idx = dp_i;
                    maximum = first[dp_i] + dotprod;
                }
            }
        }

        //void flandmark_argmax(double* smax, FLANDMARK_Options* options, const int* mapTable, FLANDMARK_PSI_SPARSE* Psi_sparse, double** q, double** g)
        //{
        void argmax(double[] smax, ref FLANDMARK_Options options, int[] mapTable, FLANDMARK_PSI_SPARSE[] Psi_sparse, List<double[]> q, List<double[]> g)
        {
            //    uint8_t M = options->M;
            byte M = options.M;

            //    // compute argmax
            //    int* indices = (int*)malloc(M * sizeof(int));
            int[] indices = new int[M];
            //    int tsize = mapTable[INDEX(1, 3, M)] - mapTable[INDEX(1, 2, M)] + 1;
            int tsize = mapTable[DefINDEX(1, 3, M)] - mapTable[DefINDEX(1, 2, M)] + 1;

            //    // left branch - store maximum and index of s5 for all positions of s1
            //    int q1_length = Psi_sparse[1].PSI_COLS;
            int q1_length = (int)Psi_sparse[1].PSI_COLS;

            //    double* s1 = (double*)calloc(2 * q1_length, sizeof(double));
            //    double* s1_maxs = (double*)calloc(q1_length, sizeof(double));
            double[] s1 = new double[2 * q1_length];
            double[] s1_maxs = new double[q1_length];

            //    for (int i = 0; i < q1_length; ++i)
            //    {
            //        // dot product <g_5, PsiGS1>
            //        flandmark_maximize_gdotprod(
            //                //s2_maxs, s2_idxs,
            //                &s1[INDEX(0, i, 2)], (double*)&s1[INDEX(1, i, 2)],
            //                q[5], g[4], options->PsiGS1[INDEX(i, 0, options->PSIG_ROWS[1])].disp,
            //                options->PsiGS1[INDEX(i, 0, options->PSIG_ROWS[1])].COLS, tsize
            //                );
            //        s1[INDEX(0, i, 2)] += q[1][i];
            //    }
            for (int i = 0; i < q1_length; ++i)
            {
                // dot product <g_5, PsiGS1>
                maximize_gdotprod(
                        //s2_maxs, s2_idxs,
                        ref s1[DefINDEX(0, i, 2)], ref s1[DefINDEX(1, i, 2)],
                        q[5], g[4], options.PsiGS1[DefINDEX(i, 0, options.PSIG_ROWS[1])].disp,
                        options.PsiGS1[DefINDEX(i, 0, options.PSIG_ROWS[1])].COLS, tsize
                        );
                s1[DefINDEX(0, i, 2)] += q[1][i];
            }

            //    for (int i = 0; i < q1_length; ++i)
            //    {
            //        s1_maxs[i] = s1[INDEX(0, i, 2)];
            //    }
            for (int i = 0; i < q1_length; ++i)
            {
                s1_maxs[i] = s1[DefINDEX(0, i, 2)];
            }

            //    // right branch (s2->s6) - store maximum and index of s6 for all positions of s2
            //    int q2_length = Psi_sparse[2].PSI_COLS;
            int q2_length = (int)Psi_sparse[2].PSI_COLS;

            //    double* s2 = (double*)calloc(2 * q2_length, sizeof(double));
            //    double* s2_maxs = (double*)calloc(q2_length, sizeof(double));
            double[] s2 = new double[2 * q2_length];
            double[] s2_maxs = new double[q2_length];

            //    for (int i = 0; i < q2_length; ++i)
            //    {
            //        // dot product <g_6, PsiGS2>
            //        flandmark_maximize_gdotprod(
            //                //s2_maxs, s2_idxs,
            //                &s2[INDEX(0, i, 2)], (double*)&s2[INDEX(1, i, 2)],
            //                q[6], g[5], options->PsiGS2[INDEX(i, 0, options->PSIG_ROWS[2])].disp,
            //                options->PsiGS2[INDEX(i, 0, options->PSIG_ROWS[2])].COLS, tsize);
            //        s2[INDEX(0, i, 2)] += q[2][i];
            //    }
            for (int i = 0; i < q2_length; ++i)
            {
                // dot product <g_6, PsiGS2>
                maximize_gdotprod(
                        //s2_maxs, s2_idxs,
                        ref s2[DefINDEX(0, i, 2)], ref s2[DefINDEX(1, i, 2)],
                        q[6], g[5], options.PsiGS2[DefINDEX(i, 0, options.PSIG_ROWS[2])].disp,
                        options.PsiGS2[DefINDEX(i, 0, options.PSIG_ROWS[2])].COLS, tsize);
                s2[DefINDEX(0, i, 2)] += q[2][i];
            }

            //    for (int i = 0; i < q2_length; ++i)
            //    {
            //        s2_maxs[i] = s2[INDEX(0, i, 2)];
            //    }
            for (int i = 0; i < q2_length; ++i)
            {
                s2_maxs[i] = s2[DefINDEX(0, i, 2)];
            }

            //    // the root s0 and its connections
            //    int q0_length = Psi_sparse[0].PSI_COLS;
            int q0_length = (int)Psi_sparse[0].PSI_COLS;

            //    double maxs0 = -FLT_MAX; int maxs0_idx = -1;
            //    double maxq10 = -FLT_MAX, maxq20 = -FLT_MAX, maxq30 = -FLT_MAX, maxq40 = -FLT_MAX, maxq70 = -FLT_MAX;
            double maxs0 = -float.MaxValue;
            int maxs0_idx = -1;
            double maxq10 = -float.MaxValue, maxq20 = -float.MaxValue, maxq30 = -float.MaxValue, maxq40 = -float.MaxValue, maxq70 = -float.MaxValue;

            //    double* s0 = (double*)calloc(M * q0_length, sizeof(double));
            double[] s0 = new double[M * q0_length];

            //    for (int i = 0; i < q0_length; ++i)
            //    {
            //        // q10
            //        maxq10 = -FLT_MAX;
            //        flandmark_maximize_gdotprod(
            //                &maxq10, &s0[INDEX(1, i, M)],
            //                s1_maxs, g[0], options->PsiGS0[INDEX(i, 0, options->PSIG_ROWS[0])].disp,
            //                options->PsiGS0[INDEX(i, 0, options->PSIG_ROWS[0])].COLS, tsize);
            //        s0[INDEX(5, i, M)] = s1[INDEX(1, (int)s0[INDEX(1, i, M)], 2)];
            //        // q20
            //        maxq20 = -FLT_MAX;
            //        flandmark_maximize_gdotprod(
            //                &maxq20, &s0[INDEX(2, i, M)],
            //                s2_maxs, g[1], options->PsiGS0[INDEX(i, 1, options->PSIG_ROWS[0])].disp,
            //                options->PsiGS0[INDEX(i, 1, options->PSIG_ROWS[0])].COLS, tsize);
            //        s0[INDEX(6, i, M)] = s2[INDEX(1, (int)s0[INDEX(2, i, M)], 2)];
            //        // q30
            //        maxq30 = -FLT_MAX;
            //        flandmark_maximize_gdotprod(
            //                &maxq30, &s0[INDEX(3, i, M)],
            //                q[3], g[2], options->PsiGS0[INDEX(i, 2, options->PSIG_ROWS[0])].disp,
            //                options->PsiGS0[INDEX(i, 2, options->PSIG_ROWS[0])].COLS, tsize);
            //        // q40
            //        maxq40 = -FLT_MAX;
            //        flandmark_maximize_gdotprod(
            //                &maxq40, &s0[INDEX(4, i, M)],
            //                q[4], g[3], options->PsiGS0[INDEX(i, 3, options->PSIG_ROWS[0])].disp,
            //                options->PsiGS0[INDEX(i, 3, options->PSIG_ROWS[0])].COLS, tsize);
            //        // q70
            //        maxq70 = -FLT_MAX;
            //        flandmark_maximize_gdotprod(
            //                &maxq70, &s0[INDEX(7, i, M)],
            //                q[7], g[6], options->PsiGS0[INDEX(i, 4, options->PSIG_ROWS[0])].disp,
            //                options->PsiGS0[INDEX(i, 4, options->PSIG_ROWS[0])].COLS, tsize);
            //        // sum q10+q20+q30+q40+q70
            //        if (maxs0 < maxq10 + maxq20 + maxq30 + maxq40 + maxq70 + q[0][i])
            //        {
            //            maxs0_idx = i;
            //            s0[INDEX(0, i, M)] = i;
            //            maxs0 = maxq10 + maxq20 + maxq30 + maxq40 + maxq70 + q[0][i];
            //        }
            //    }
            for (int i = 0; i < q0_length; ++i)
            {
                // q10
                maxq10 = -float.MaxValue;
                maximize_gdotprod(
                        ref maxq10, ref s0[DefINDEX(1, i, M)],
                        s1_maxs, g[0], options.PsiGS0[DefINDEX(i, 0, options.PSIG_ROWS[0])].disp,
                        options.PsiGS0[DefINDEX(i, 0, options.PSIG_ROWS[0])].COLS, tsize);
                s0[DefINDEX(5, i, M)] = s1[DefINDEX(1, (int)s0[DefINDEX(1, i, M)], 2)];
                // q20
                maxq20 = -float.MaxValue;
                maximize_gdotprod(
                        ref maxq20, ref s0[DefINDEX(2, i, M)],
                        s2_maxs, g[1], options.PsiGS0[DefINDEX(i, 1, options.PSIG_ROWS[0])].disp,
                        options.PsiGS0[DefINDEX(i, 1, options.PSIG_ROWS[0])].COLS, tsize);
                s0[DefINDEX(6, i, M)] = s2[DefINDEX(1, (int)s0[DefINDEX(2, i, M)], 2)];
                // q30
                maxq30 = -float.MaxValue;
                maximize_gdotprod(
                        ref maxq30, ref s0[DefINDEX(3, i, M)],
                        q[3], g[2], options.PsiGS0[DefINDEX(i, 2, options.PSIG_ROWS[0])].disp,
                        options.PsiGS0[DefINDEX(i, 2, options.PSIG_ROWS[0])].COLS, tsize);
                // q40
                maxq40 = -float.MaxValue;
                maximize_gdotprod(
                        ref maxq40, ref s0[DefINDEX(4, i, M)],
                        q[4], g[3], options.PsiGS0[DefINDEX(i, 3, options.PSIG_ROWS[0])].disp,
                        options.PsiGS0[DefINDEX(i, 3, options.PSIG_ROWS[0])].COLS, tsize);
                // q70
                maxq70 = -float.MaxValue;
                maximize_gdotprod(
                        ref maxq70, ref s0[DefINDEX(7, i, M)],
                        q[7], g[6], options.PsiGS0[DefINDEX(i, 4, options.PSIG_ROWS[0])].disp,
                        options.PsiGS0[DefINDEX(i, 4, options.PSIG_ROWS[0])].COLS, tsize);
                // sum q10+q20+q30+q40+q70
                if (maxs0 < maxq10 + maxq20 + maxq30 + maxq40 + maxq70 + q[0][i])
                {
                    maxs0_idx = i;
                    s0[DefINDEX(0, i, M)] = i;
                    maxs0 = maxq10 + maxq20 + maxq30 + maxq40 + maxq70 + q[0][i];
                }
            }

            //    // get indices
            //    for (int i = 0; i < M; ++i)
            //    {
            //        indices[i] = (int)s0[INDEX(0, maxs0_idx, M) + i] + 1;
            //    }
            for (int i = 0; i < M; ++i)
            {
                indices[i] = (int)s0[DefINDEX(0, maxs0_idx, M) + i] + 1;
            }

            //    // cleanup temp variables
            //    free(s0);
            //    free(s1); free(s1_maxs);
            //    free(s2); free(s2_maxs);

            //    // convert 1D indices to 2D coordinates of estimated positions
            //    //int * optionsS = &options->S[0];
            //    const int* optionsS = options->S;
            //    for (int i = 0; i < M; ++i)
            //    {
            //        int rows = optionsS[INDEX(3, i, 4)] - optionsS[INDEX(1, i, 4)] + 1;
            //        smax[INDEX(0, i, 2)] = float(COL(indices[i], rows) + optionsS[INDEX(0, i, 4)]);
            //        smax[INDEX(1, i, 2)] = float(ROW(indices[i], rows) + optionsS[INDEX(1, i, 4)]);
            //    }
            for (int i = 0; i < M; ++i)
            {
                int rows = options.S[DefINDEX(3, i, 4)] - options.S[DefINDEX(1, i, 4)] + 1;
                smax[DefINDEX(0, i, 2)] = (float)(DefCOL(indices[i], rows) + options.S[DefINDEX(0, i, 4)]);
                smax[DefINDEX(1, i, 2)] = (float)(DefROW(indices[i], rows) + options.S[DefINDEX(1, i, 4)]);
            }

            //    free(indices);
            //}
        }
    }
}
