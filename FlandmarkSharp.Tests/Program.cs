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
            var modeldata = new LibFlandmark(Path.Combine(Environment.CurrentDirectory, @"data\flandmark_model.dat"));

            ImageTest(modeldata);
        }

        static void VideoTest(LibFlandmark modeldata)
        {
            using(CascadeClassifier cascade = new CascadeClassifier("haarcascade_frontalface_alt.xml"))
            using (VideoCapture capture = new VideoCapture())
            {
                capture.Open(1);
                while (true)
                {
                    using (Mat frame = new Mat())
                    {
                        if (capture.Read(frame) && !frame.Empty())
                        {
                            using (Mat m = frame.Clone())
                            {
                                m.Resize(new Size(m.Width * 0.2, m.Height * 0.2));

                                var result = cascade.DetectMultiScale(m);
                                if (result != null && result.Length > 0)
                                {
                                    Rect r = result[0];
                                    Cv2.Rectangle(frame, r, Scalar.Green);

                                    DetectAndDraw(modeldata, frame, new int[] { Math.Max(r.X, 0), Math.Max(0, r.Y), Math.Min(r.X + r.Width, m.Width - 1), Math.Min(r.Y + r.Height, m.Height - 1) });
                                }
                            }

                            Cv2.ImShow("cam", frame);
                            Cv2.WaitKey(10);
                        }
                    }
                }
            }
        }

        static void ImageTest(LibFlandmark modeldata)
        {
            string dir = Path.Combine(Environment.CurrentDirectory, @"data\");
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

                    DetectAndDraw(modeldata, m, bb);

                    Cv2.ImShow("img", m);
                    char c = (char)Cv2.WaitKey(3000);
                    if (c == 'e')
                        return;
                }
            }
        }

        static void DetectAndDraw(LibFlandmark f, Mat m, int[] box)
        {
            Mat gray = m.CvtColor(ColorConversionCodes.BGR2GRAY);

            int[] bb = box;

            var result = f.Detect(gray, bb);

            if (result != null)
            {
                for (int i = 0; i < result.Length / 2; i++)
                {
                    Scalar color = Scalar.Red;
                    if (i == 0)
                        color = Scalar.Magenta;
                    Cv2.Circle(m, new Point(result[i * 2], result[i * 2 + 1]), 2, color, 4);
                }
            }

            Cv2.Rectangle(m, new Rect(bb[0], bb[1], bb[2] - bb[0], bb[3] - bb[1]), Scalar.Blue, 2);
        }
    }
}
