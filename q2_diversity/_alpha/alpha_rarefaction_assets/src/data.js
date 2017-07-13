export default function setupData(data) {
  const [xAxisLabel, yAxisLabel] = data.data.columns;

  let minX = Infinity;
  let maxX = 0;
  let minY = Infinity;
  let maxY = 0;

  data.data.data.forEach((d) => {
    const [x, y] = d;
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    if (y < minY) minY = y;
    if (y > maxY) maxY = y;
  });

  return {
    data: data.data.data,
    xAxisLabel,
    yAxisLabel,
    minX,
    maxX,
    minY,
    maxY,
  };
}
