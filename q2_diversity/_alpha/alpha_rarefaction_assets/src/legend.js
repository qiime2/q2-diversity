import { select } from 'd3';
import { curData, toggle, toggleAll } from './data';

export default function appendLegendKey(legend, entry, ly, c, color) {
  // line toggle in the legend
  console.log(entry, ':', curData[entry]);
  const rect = legend.append('rect')
    .attr('id', `idrect${entry}`)
    .attr('class', 'legend rect')
    .attr('x', 0)
    .attr('y', ly - 2.5)
    .attr('width', 15)
    .attr('height', 5)
    .attr('stroke', 'darkGrey')
    .attr('fill', () => curData[entry].line);
  rect.on('click', () => {
    const newColor = curData[entry].line === c ? 'white' : c;
    if (entry === 'Select%20All') {
      toggleAll(false, true, color);
      for (const key of curData) {
        const thatRect = select(`#idrect${key}`);
        thatRect.attr('fill', curData[key].line);
      }
    } else {
      toggle(entry, null, newColor);
      rect.attr('fill', curData[entry].line);
    }
  });
  // dot toggle in the legend
  const dot = legend.append('circle')
    .attr('id', `iddot${entry}`)
    .attr('class', 'legend circle')
    .attr('cx', 30)
    .attr('cy', ly)
    .attr('r', 5)
    .attr('stroke', 'darkGrey')
    .attr('fill', () => curData[entry].dots);
  dot.on('click', () => {
    const newColor = curData[entry].dots === c ? 'white' : c;
    if (entry === 'Select%20All') {
      toggleAll(true, false, color);
      for (const key of curData) {
        const thatCircle = select(`#iddot${key}`);
        thatCircle.attr('fill', curData[key].dots);
      }
    } else {
      toggle(entry, newColor, null);
      dot.attr('fill', () => curData[entry].dots);
    }
  });
  // text for key in the legend
  legend.append('text')
      .attr('class', 'legend')
      .attr('x', 40)
      .attr('y', ly + 5)
      .attr('font', '10px sans-serif')
      .text(decodeURI(entry));
}

