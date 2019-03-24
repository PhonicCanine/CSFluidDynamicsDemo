using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Threading;

namespace CSFluidDynamics
{
    public class simPoint
    {
        public decimal x { get; set; }
        public decimal y { get; set; }
        public static vector operator +(simPoint a, simPoint b)
        {
            return new vector() { dx = a.x + b.x, dy = a.y + b.y };
        }
        public static simPoint operator +(simPoint a, vector b)
        {
            return new simPoint() { x = a.x + b.dx, y = a.y + b.dy };
        }

        public static vector operator -(simPoint a, simPoint b)
        {
            return new vector() { dx = a.x - b.x, dy = a.y - b.y };
        }
        public static simPoint operator -(simPoint a, vector b)
        {
            return new simPoint() { x = a.x - b.dx, y = a.y - b.dy };
        }

        public vector asVector {
            get
            {
                return new vector(x, y);
            } }

        public static simPoint fromPointCollection(IEnumerable<simPoint> points)
        {
            vector zero = new vector();
            int count = 0;
            foreach (simPoint v in points)
            {
                zero += v.asVector;
                count += 1;
            }
            return (zero / count).asPoint;
        }

        public simPoint()
        {

        }

        public simPoint(decimal posx, decimal posy)
        {
            x = posx;
            y = posy;
        }
    }

    public class vector
    {
        public vector()
        {

        }
        public vector(decimal x, decimal y)
        {
            dx = x;
            dy = y;
        }
        public vector(double magnitude, double angle)
        {
            double y = Math.Sin(angle) * magnitude;
            double x = Math.Cos(angle) * magnitude;
            dx = (decimal)x;
            dy = (decimal)y;
        }

        public decimal dx { get; set; }
        public decimal dy { get; set; }
        public static vector operator *(vector a, decimal b)
        {
            return new vector() { dx = a.dx * b, dy = a.dy * b };
        }
        public static vector operator /(vector a, decimal b)
        {
            if (b == 0)
                return new vector();
            return new vector() { dx = a.dx / b, dy = a.dy / b };
        }
        public static vector operator +(vector a, vector b)
        {
            return new vector() { dx = a.dx + b.dx, dy = a.dy + b.dy };
        }
        public static vector operator -(vector a, vector b)
        {
            return new vector() { dx = a.dx - b.dx, dy = a.dy - b.dy };
        }

        //https://stackoverflow.com/questions/4124189/performing-math-operations-on-decimal-datatype-in-c
        // x - a number, from which we need to calculate the square root
        // epsilon - an accuracy of calculation of the root from our number.
        // The result of the calculations will differ from an actual value
        // of the root on less than epslion.
        public static decimal Sqrt(decimal x, decimal epsilon = 0.0M)
        {
            if (x < 0) throw new OverflowException("Cannot calculate square root from a negative number");

            decimal current = (decimal)Math.Sqrt((double)x), previous;
            do
            {
                previous = current;
                if (previous == 0.0M) return 0;
                current = (previous + x / previous) / 2;
            }
            while (Math.Abs(previous - current) > epsilon);
            return current;
        }

        public decimal magnitude { get
            {
                try
                {
                    return Sqrt((dx * dx) + (dy * dy));
                }
                catch (System.OverflowException)
                {
                    return 0;
                }
                
            } }
        public vector normalised { get
            {
                try
                {
                    vector vec = new vector() { dx = this.dx, dy = this.dy };
                    vec = vec / vec.magnitude;
                    return vec;
                }
                catch (System.DivideByZeroException)
                {
                    return vector.zero;
                }
                
            } }
        public double angle { get
            {
                return Math.Tanh((double)dy / (double)dx);
            } }

        public static vector zero
        {
            get
            {
                return new vector();
            }
        }

        public simPoint asPoint
        {
            get
            {
                return new simPoint() { x = dx, y = dy };
            }
        }

        public static vector fromVectorCollection(IEnumerable<vector> vectors)
        {
            vector zero = new vector();
            int count = 0;
            foreach (vector v in vectors)
            {
                zero += v;
                count += 1;
            }
            return zero / count;
        }
    }

    public class particle
    {
        public decimal repellanceCoefficient { get; set; }
        public simPoint position { get; set; }
        public vector momentum { get; set; }
        public vector momentumChange { get; set; }

        private Guid id;
        public Guid identifier { get
            {
                return id;
            }
        }

        //some averagey number
        public decimal charge { get
            {
                return repellanceCoefficient;
            } }

        public particle()
        {
            repellanceCoefficient = -10;
            id = Guid.NewGuid();
        }

        public particle copy
        {
            get
            {
                return new particle() { momentum=momentum,repellanceCoefficient=repellanceCoefficient,position=position,id=id };
            }
        }

        public particle copyWithNewID
        {
            get
            {
                return new particle() { momentum = momentum, repellanceCoefficient = repellanceCoefficient, position = position, id = Guid.NewGuid() };
            }
        }
    }

    public class simulationCell
    {
        #region surroundingCells
        public simulationCell topLeft
        {
            get
            {
                return simulationCells[0][0];
            }
            set
            {
                simulationCells[0][0] = value;
            }
        }

        public simulationCell top
        {
            get
            {
                return simulationCells[1][0];
            }
            set
            {
                simulationCells[1][0] = value;
            }
        }

        public simulationCell topRight
        {
            get
            {
                return simulationCells[2][0];
            }
            set
            {
                simulationCells[2][0] = value;
            }
        }

        public simulationCell right
        {
            get
            {
                return simulationCells[2][1];
            }
            set
            {
                simulationCells[2][1] = value;
            }
        }

        public simulationCell bottomRight
        {
            get
            {
                return simulationCells[2][2];
            }
            set
            {
                simulationCells[2][2] = value;
            }
        }

        public simulationCell bottom
        {
            get
            {
                return simulationCells[1][2];
            }
            set
            {
                simulationCells[1][2] = value;
            }
        }

        public simulationCell bottomLeft
        {
            get
            {
                return simulationCells[0][2];
            }
            set
            {
                simulationCells[0][2] = value;
            }
        }

        public simulationCell left
        {
            get
            {
                return simulationCells[0][1];
            }
            set
            {
                simulationCells[0][1] = value;
            }
        }

        #endregion

        //indexed x,y
        //x + -> left to right
        //y + -> top to bottom
        public simulationCell[][] simulationCells = new simulationCell[3][] { new simulationCell[3], new simulationCell[3], new simulationCell[3] };

        public simulationCell this[int x, int y]
        {
            get
            {
                if (x == y && x == 0)
                {
                    return this;
                }
                else
                {
                    return simulationCells[x + 1][y + 1];
                }
            }
        }
        public simulationCell this[int x]
        {
            get
            {
                int col = x % 3;
                int row = (int)(x / 3);
                return this[row - 1, col - 1];
            }
        }

        private List<particle> containedParticles { get; set; }
        private List<simPoint> containedGeometry { get; set; }

        public simPoint topLeftBoundary { get; set; }
        public simPoint bottomRightBoundary { get; set; }

        public List<particle> getParticles { get
            {
                return containedParticles;
            } }

        public bool isPointContained(simPoint pt)
        {
            if (pt.x < bottomRightBoundary.x && pt.x > topLeftBoundary.x && pt.y < bottomRightBoundary.y && pt.y > topLeftBoundary.y)
                return true;
            return false;
        }

        private List<particle> particlesToAddToList;
        private List<particle> particlesToRemoveFromList;
        public void addParticlesToAddList(IEnumerable<particle> ps)
        {
            if (particlesToAddToList == null)
                particlesToAddToList = new List<particle>();
            particlesToAddToList.AddRange(ps);
        }

        //later use this to provide simple repulsion of a particle from a chunk
        public decimal chunkDensity { get
            {
                if (!activeCell)
                    return defaultDensity;
                return (decimal)containedParticles.Count / ((bottomRightBoundary.x - topLeftBoundary.x) * (bottomRightBoundary.y - topLeftBoundary.y));
            } }

        private bool activeCell = true;
        private decimal defaultDensity = 0;

        public bool active
        {
            get
            {
                return activeCell;
            }
        }

        public simulationCell(int defaultDensity_, simPoint tlBoundary, simPoint brBoundary)
        {
            activeCell = false;
            defaultDensity = defaultDensity_;
            topLeftBoundary = tlBoundary;
            bottomRightBoundary = brBoundary;
        }

        private particle initialPrototype;
        private vector initialVectorPrototype;

        public simulationCell(simPoint tlBoundary, simPoint brBoundary, List<simPoint> geometry, int numberOfParticles, particle prototype, vector startingVector = null)
        {
            topLeftBoundary = tlBoundary;
            bottomRightBoundary = brBoundary;
            containedGeometry = geometry;

            initialPrototype = prototype;
            initialVectorPrototype = startingVector;

            spawnParticles(numberOfParticles,prototype,startingVector);
        }

        public void addNewParticles(int num)
        {
            if (activeCell)
                spawnParticles(num, initialPrototype, initialVectorPrototype);
        }

        private void spawnParticles(int particles, particle proto, vector v)
        {
            var height = bottomRightBoundary.y - topLeftBoundary.y;
            var width = bottomRightBoundary.x - topLeftBoundary.x;
            if (height == 0 || width == 0)
                throw new Exception("Can't define simulation zone with zero area!");
            decimal ratio = height / width;
            decimal inverseRatio = 1 / ratio;
            decimal numberToMultiply = vector.Sqrt(particles);

            int wide = (int)(inverseRatio * numberToMultiply);
            int high = (int)(ratio * numberToMultiply);

            if (containedParticles == null)
                containedParticles = new List<particle>();

            for (int x = 0; x < wide; x++)
            {
                for (int y = 0; y < high; y++)
                {
                    particle p = proto.copyWithNewID;
                    decimal posx = topLeftBoundary.x + (width / wide) * x;
                    decimal posy = topLeftBoundary.y + (height / high) * y;
                    p.position = new simPoint(posx, posy);
                    if (v != null)
                        p.momentum = v;

                    containedParticles.Add(p);
                }
            }
        }

        private List<particle> particlesToAdd;
        private List<Guid> particlesToRemove;

        //o(n^2) -- terrible efficiency overall but it might work
        public void runSimulation()
        {

            if (!this.active)
                return;
            lock (this)
            {
                if (particlesToAdd != null)
                    containedParticles.AddRange(particlesToAdd);
                lock (containedParticles)
                {
                    List<particle> toRemoveNow = new List<particle>();
                    if (particlesToRemove != null)
                        foreach (Guid pid in particlesToRemove)
                            foreach (particle p in containedParticles)
                                if (p.identifier == pid)
                                    toRemoveNow.Add(p);
                    foreach (particle p in toRemoveNow)
                        containedParticles.Remove(p);
                }
                

                particlesToRemove = new List<Guid>();
                particlesToAdd = new List<particle>();
            }

            simPoint av = simPoint.fromPointCollection(new[] { topLeftBoundary, bottomRightBoundary });

            List<particle> ps = new List<particle>();
            foreach (particle p in containedParticles)
            {
                if (p != null)
                    ps.Add(p);
            }
                //p == null ? (() => { }). : ps.Add(p);
            containedParticles = ps;

            for (int i = 0; i < containedParticles.Count; i++)
            {
                vector momentumToAdd = new vector();
                for (int z = 0; z < 9; z++)
                {
                    if (this[z] != null && this[z].active)
                    {
                        for (int j = 0; j < this[z].containedParticles.Count; j++)
                        {
                            if ((i != j || z != 4) && (j < this[z].containedParticles.Count && i < containedParticles.Count))
                            {
                                try
                                {
                                    var vecDiff = containedParticles[i].position - this[z].containedParticles[j].position;
                                    vecDiff = vecDiff.normalised * simulation.forceBetweenParticles(this[z].containedParticles[j], containedParticles[i]);
                                    momentumToAdd += vecDiff;
                                }
                                catch (ArgumentOutOfRangeException)
                                {
                                    //don't know how this can go wrong, but, it can
                                }
                                catch (NullReferenceException)
                                {
                                    //Also don't know how this happens, but sometimes it does
                                }
                                
                            }
                        }
                    }
                    
                }
                
                for (int gp = 0; gp < containedGeometry.Count; gp++)
                {
                    var vecDiff = containedParticles[i].position - containedGeometry[gp];
                    vecDiff = vecDiff.normalised * simulation.forceBetweenParticles(new particle() { position = containedGeometry[gp], momentum = new vector(), repellanceCoefficient = simulation.solidCoefficient}, containedParticles[i]);
                    momentumToAdd += vecDiff;
                }

                simPoint normalised = (containedParticles[i].position - av).asPoint;

                //needs to push back harder if the density is lower in the current chunk than the other chunk. Should probably use distance for this? or charge? idk
                //trying with charge first because linear relations are more easy to understand, but distance is going to be the one that needs to change, because density is proportional to distance, not charge.
                //not going to work, because there's no way of knowing where a particle lies in relation to the simulationchunk.
                decimal density = chunkDensity;

                //top
                simulationCell top = this[0, -1];
                if (top != null)
                {
                    decimal topDensity = top.chunkDensity;
                    var downVec = new vector(0, 1M);
                    if (!top.active)
                    {
                        momentumToAdd += (downVec * containedParticles[i].momentum.magnitude) / (normalised).asVector.magnitude;
                    }
                }
                

                //right
                simulationCell right = this[1, 0];
                if (right != null)
                {
                    decimal rightDensity = right.chunkDensity;
                    var leftVec = new vector(-1M, 0);
                    if (!right.active)
                    {
                        momentumToAdd += (leftVec * containedParticles[i].momentum.magnitude) / (normalised).asVector.magnitude;
                    }
                }
                

                //bottom
                simulationCell bottom = this[0, 1];
                if (bottom != null)
                {
                    decimal bottomDensity = bottom.chunkDensity;
                    var upVec = new vector(0, -1M);
                    if (!bottom.active)
                    {
                        momentumToAdd += (upVec * containedParticles[i].momentum.magnitude) / (normalised).asVector.magnitude;
                    }
                }
                

                //left
                simulationCell left = this[-1, 0];
                if (left != null)
                {
                    decimal leftDensity = left.chunkDensity;
                    var rightVec = new vector(1M, 0);
                    if (!left.active)
                    {
                        //momentumToAdd += (rightVec * containedParticles[i].momentum.magnitude) / (normalised).asVector.magnitude;
                    }
                }

                containedParticles[i].momentumChange = momentumToAdd;

                containedParticles[i].momentum += momentumToAdd;
            }

            for (int i = 0; i < containedParticles.Count; i++)
            {
                containedParticles[i].position += containedParticles[i].momentum * simulation.simulationDelta;
                if (this.isPointContained(containedParticles[i].position) == false)
                {
                    bool relocated = false;

                    for (int z = 0; z < 9; z++)
                    {
                        if (this[z] != null)
                        {
                            if (this[z].isPointContained(containedParticles[i].position) && z != 4)
                            {
                                this[z].addParticlesToAddList(new[] { containedParticles[i] });
                                if (particlesToRemoveFromList == null)
                                    particlesToRemoveFromList = new List<particle>();
                                particlesToRemoveFromList.Add(containedParticles[i]);
                                relocated = true;
                            }
                        }

                    }


                    if (!relocated)
                    {
                        if (particlesToRemoveFromList == null)
                            particlesToRemoveFromList = new List<particle>();
                        particlesToRemoveFromList.Add(containedParticles[i]);
                    }
                }
                
            }

        }

        public void cleanup()
        {

            lock (this)
            {
                particlesToAdd = particlesToAddToList;
                particlesToAddToList = null;
                particlesToRemove = new List<Guid>();
                if (particlesToRemoveFromList != null)
                    foreach (particle p in particlesToRemoveFromList)
                        particlesToRemove.Add(p.identifier);
                particlesToRemoveFromList = null;
            }
        }
    }

    //needs to start with momentum in one direction or another
    //needs to split simulation into chunks most likely, because then search should be easier
    public class simulation
    {
        List<simulationCell> cells;

        public const decimal k = 0.9M;
        public const decimal q = 2;
        public const decimal scale = 1M / 1M;
        public const decimal simulationDelta = 1M / 50M;

        public const decimal solidCoefficient = -1;
        public const decimal airCoefficient = -6;

        public static decimal forceBetweenParticles(particle acting, particle actedUpon)
        {
            var posDiff = actedUpon.position - acting.position;
            decimal dist = posDiff.magnitude * scale;
            if (dist == 0)
                dist = 0.000000001M;
            return (k * acting.charge * q * actedUpon.charge * q) / (dist * dist);

        }

        private List<simulationCell> spawnCells;

        private int spawncount = 0;

        //air will always flow right to left
        public simulation(List<simPoint> geometryPoints, int width = 700, int height = 300, int cellsHigh = 6, int cellsWide = 16, decimal baseDensity = 0.0125M, vector startingVelocity = null, particle p = null)
        {

            spawnCells = new List<simulationCell>();

            int realWidth = width * 2;
            int realCellsWide = cellsWide + 2;//* 2;
            if (startingVelocity == null)
                startingVelocity = new vector(-0.1M, 0);
            if (p == null)
                p = new particle() { momentum = startingVelocity, repellanceCoefficient = airCoefficient, position = vector.zero.asPoint };

            if (cells == null)
                cells = new List<simulationCell>();


            decimal cellWidth = (decimal)(width / cellsWide);
            decimal cellHeight = (decimal)(height / cellsHigh);

            int particlesPerCell = (int)(cellWidth * cellHeight * baseDensity);
            spawncount = particlesPerCell;

            Dictionary<Tuple<int, int>, simulationCell> cellDict = new Dictionary<Tuple<int, int>, simulationCell>();

            for (int x = -1; x <= realCellsWide; x++)
            {
                for (int y = -1; y <= cellsHigh; y++)
                {
                    simPoint topLeftB = new simPoint(x * cellWidth, y * cellHeight);
                    simPoint bottomRightB = new simPoint((x + 1) * cellWidth, (y + 1) * cellHeight);

                    //make real simulation cells
                    simulationCell c;
                    if (y >= 0 && y < cellsHigh)
                    {
                        if (x > cellsWide && x < realCellsWide)
                        {
                            c = new simulationCell(topLeftB, bottomRightB, new List<simPoint>(), particlesPerCell, p);
                            spawnCells.Add(c);
                        }
                        else if (x >= 0)
                        {
                            List<simPoint> ptsRelevant = geometryPoints.FindAll((s) => {
                                if (s.x <= bottomRightB.x && s.x >= topLeftB.x && s.y >= topLeftB.y && s.y < bottomRightB.y)
                                {
                                    return true;
                                }
                                else
                                {
                                    return false;
                                }
                            });
                            c = new simulationCell(topLeftB, bottomRightB, ptsRelevant, 0, p);
                            
                        }
                        else
                        {
                            c = new simulationCell(particlesPerCell / 2, topLeftB, bottomRightB);
                        }
                    }
                    else
                    {
                        c = new simulationCell(particlesPerCell / 2, topLeftB, bottomRightB);
                    }
                    cellDict[new Tuple<int, int>(x, y)] = c;
                    cells.Add(c);
                    
                }
            }

            for (int x = 0; x < realCellsWide; x++)
            {
                for (int y = 0; y < cellsHigh; y++)
                {
                    cellDict[new Tuple<int, int>(x, y)].bottom = cellDict[new Tuple<int, int>(x, y + 1)];
                    cellDict[new Tuple<int, int>(x, y)].top = cellDict[new Tuple<int, int>(x, y - 1)];
                    cellDict[new Tuple<int, int>(x, y)].right = cellDict[new Tuple<int, int>(x + 1, y)];
                    cellDict[new Tuple<int, int>(x, y)].left = cellDict[new Tuple<int, int>(x - 1, y)];
                    cellDict[new Tuple<int, int>(x, y)].bottomRight = cellDict[new Tuple<int, int>(x + 1, y + 1)];
                    cellDict[new Tuple<int, int>(x, y)].bottomLeft = cellDict[new Tuple<int, int>(x - 1, y + 1)];
                    cellDict[new Tuple<int, int>(x, y)].topRight = cellDict[new Tuple<int, int>(x + 1, y - 1)];
                    cellDict[new Tuple<int, int>(x, y)].topLeft = cellDict[new Tuple<int, int>(x - 1, y - 1)];
                }
            }

        }

        //loop through all simulationcells, calling runSimulation once, then cleanup
        //use parallel.for
        public void simulationCycle()
        {

            for (int x = 0; x < spawnCells.Count; x++) {
                if (spawnCells[x].getParticles.Count < spawncount / 2)
                    spawnCells[x].addNewParticles(spawncount);
            }

            Parallel.For(0, cells.Count, (int x) =>
            {
                cells[x].runSimulation();
                
                //System.Diagnostics.Debug.WriteLine("Did cell " + x.ToString());
            });

            for (int x = 0; x < cells.Count; x++)
            {
                cells[x].cleanup();
                //System.Diagnostics.Debug.WriteLine("Cleaned cell " + x.ToString());
            }
        }

        public List<particle> fetchAllParticles()
        {
            List<particle> ptc = new List<particle>();
            foreach (simulationCell cell in cells)
            {
                var particles = cell.getParticles;
                if (particles != null)
                    ptc.AddRange(cell.getParticles);
            }
            return ptc;
        }

        public List<particle> runSimulationAndFetchParticles()
        {
            simulationCycle();
            return fetchAllParticles();
        }
    }

    public class simulationVolume
    {
        public List<simPoint> objectEdges { get; set; }
    }

    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
            
        }

        public simulation simulator;

        public List<simPoint> findEdges()
        {
            var render = new RenderTargetBitmap((int)SimulationCanvas.ActualWidth, (int)SimulationCanvas.ActualHeight, 96, 96, PixelFormats.Pbgra32);
            render.Render(SimulationCanvas);

            MemoryStream stream = new MemoryStream();
            BitmapEncoder encoder = new BmpBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(render));
            encoder.Save(stream);
            var fImageData = stream.ToArray();
            Bitmap bitmap = new Bitmap(stream);

            List<simPoint> edgeParts = new List<simPoint>();
            List<simPoint> allParts = new List<simPoint>();

            for (int x = 0; x < bitmap.Width; x++)
            {
                for (int y = 0; y < bitmap.Height; y++)
                {
                    if (bitmap.GetPixel(x, y) != System.Drawing.Color.FromArgb(255,0,0,0))
                    {
                        allParts.Add(new simPoint() { x = x, y = y });

                        if (bitmap.GetPixel(x + 1, y) == System.Drawing.Color.FromArgb(255, 0, 0, 0) || bitmap.GetPixel(x - 1, y) == System.Drawing.Color.FromArgb(255, 0, 0, 0) || bitmap.GetPixel(x, y + 1) == System.Drawing.Color.FromArgb(255, 0, 0, 0) || bitmap.GetPixel(x, y - 1) == System.Drawing.Color.FromArgb(255, 0, 0, 0))
                        {
                            var item = new System.Windows.Shapes.Ellipse();
                            item.Fill = new System.Windows.Media.SolidColorBrush(System.Windows.Media.Colors.Red);
                            item.Width = 4;
                            item.Height = 4;
                            item.Margin = new Thickness(x, y, 0, 0);
                            SimulationCanvas.Children.Add(item);
                            edgeParts.Add(new simPoint() {x = x, y = y});
                        }

                        
                    }
                    else
                    {

                    }
                }
            }

            //return edgeParts;
            return allParts;
        }

        private void DoThing(object sender, RoutedEventArgs e)
        {
            var geometry = findEdges();
            simulator = new simulation(geometry);
            StepButton.IsEnabled = true;
        }

        bool canClick = true;

        private void drawSimulationState(List<particle> ps)
        {
            SimulationCanvasForeground.Children.Clear();

            decimal lowestVectorMagnitude = decimal.MaxValue;
            decimal highestVectorMagnitute = 0M;

            foreach (var p in ps)
            {
                if (p.momentumChange.magnitude < lowestVectorMagnitude)
                {
                    lowestVectorMagnitude = p.momentumChange.magnitude;
                }

                if (p.momentumChange.magnitude > highestVectorMagnitute)
                {
                    highestVectorMagnitute = p.momentumChange.magnitude;
                }
            }

            foreach (var p in ps)
            {
                var item = new System.Windows.Shapes.Ellipse();

                decimal percentile = (p.momentumChange.magnitude - lowestVectorMagnitude) / highestVectorMagnitute;

                byte red = (byte)(255 * percentile);
                byte blue = (byte)(255 * (1M - percentile));

                item.Fill = new System.Windows.Media.SolidColorBrush(System.Windows.Media.Color.FromArgb(255,red,0,blue));
                item.Width = 4;
                item.Height = 4;
                item.Margin = new Thickness((double)p.position.x, (double)p.position.y, 0, 0);
                SimulationCanvasForeground.Children.Add(item);
            }

            Title = "Particle Simulator: " + ps.Count + " particles displayed";
        }

        List<List<particle>> history = new List<List<particle>>();

        private List<particle> copy(List<particle> input)
        {
            List<particle> outp = new List<particle>();
            foreach (particle p in input)
            {
                outp.Add(new particle() { momentumChange = p.momentumChange, position = p.position });
            }
            return outp;
        }

        private void SimulationStep(object sender, RoutedEventArgs e)
        {

            if (!canClick)
                return;

            var t = new Thread(() => {
                List<particle> ps = new List<particle>();
                int max = 5000;

                for (var i = 0; i < max; i++)
                {
                    ps = simulator.runSimulationAndFetchParticles();
                    this.Dispatcher.Invoke(() =>
                    {
                        ProgressBar1.Maximum = max;
                        ProgressBar1.Minimum = 0;
                        ProgressBar1.Value = i;
                        drawSimulationState(ps);
                    });
                    
                    history.Add(copy(ps));
                }

                this.Dispatcher.Invoke(() =>
                {
                    PlaybackButton.IsEnabled = true;
                    canClick = true;
                });

                });

            canClick = false;
            PlaybackButton.IsEnabled =false;
            t.Start();

        }

        private int replayPosition = 0;

        private void Playback(object sender, RoutedEventArgs e)
        {
            replayPosition = 0;

            DispatcherTimer t = new DispatcherTimer(new TimeSpan(0, 0, 0, 0, 20), priority: DispatcherPriority.Render, (a, b) => {
                if (replayPosition < history.Count)
                {
                    drawSimulationState(history[replayPosition]);
                    replayPosition += 1;
                }
                else
                {
                    ((DispatcherTimer)a).Stop();
                    PlaybackButton.IsEnabled = true;
                }
            }, this.Dispatcher);

            PlaybackButton.IsEnabled = false;
            t.Start();

        }
    }
}
