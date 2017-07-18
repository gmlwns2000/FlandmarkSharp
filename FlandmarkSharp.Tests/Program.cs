using OpenCvSharp;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlandmarkSharp.Tests
{
    class Program
    {
        static void Main(string[] args)
        {
            var modeldata = new LibFlandmark(@"C:\Users\AinL\Documents\Visual Studio 2017\Projects\flandmark-master\data\flandmark_model.dat");
            string dir = @"C:\Users\AinL\Documents\Visual Studio 2017\Projects\flandmark-master\data\images";
            FileInfo[] files = new DirectoryInfo(dir).GetFiles();
            foreach (FileInfo fi in files)
            {
                if (fi.Name.EndsWith(".jpg"))
                {
                    Mat m = Cv2.ImRead(fi.FullName);

                    FileInfo det = new FileInfo(Path.Combine(fi.Directory.FullName, Path.GetFileNameWithoutExtension(fi.FullName) + ".det"));
                    string content = File.ReadAllText(det.FullName);
                    string[] spl = content.Split(',');
                    int[] bb = new int[4];
                    for (int i = 0; i < 4; i++)
                    {
                        bb[i] = Convert.ToInt32(spl[i]);
                    }

                    DetectAndShow(modeldata, m, bb);

                    Cv2.ImShow("img", m);
                    Cv2.WaitKey(3000);
                }
            }
        }

        static void DetectAndShow(LibFlandmark f, Mat m, int[] box)
        {
            Mat gray = m.CvtColor(ColorConversionCodes.BGR2GRAY);
            int[] bb = box;
            var result = f.Detect(gray, bb);
            for (int i = 0; i < result.Length / 2; i++)
            {
                Scalar color = Scalar.Red;
                if (i == 0)
                    color = Scalar.Purple;
                Cv2.Circle(m, new Point(result[i * 2], result[i * 2 + 1]), 2, color, 4);
            }
            Cv2.Rectangle(m, new Rect(bb[0], bb[1], bb[2] - bb[0], bb[3] - bb[1]), Scalar.Blue, 2);
        }
    }
}
