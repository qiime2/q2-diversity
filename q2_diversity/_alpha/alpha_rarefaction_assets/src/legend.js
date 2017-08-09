import { select } from 'd3';
import { curData, toggle } from './data';

export default function appendLegendKey(legend, entry, ly, color) {
  // line toggle in the legend
  const all = 'Select%20All';
  const rect = legend.append('rect')
    .attr('id', `rect${entry}`)
    .attr('class', 'legend rect')
    .attr('x', 0)
    .attr('y', ly - 2.5)
    .attr('width', 15)
    .attr('height', 5)
    .attr('stroke', 'darkGrey')
    .attr('fill', curData[entry].line);
  rect.on('click', () => {
    const becomeFull = curData[entry].lineOpacity === 0;
    if (entry === all) {
      for (const key of Object.keys(curData).slice(1)) {
        const thatRect = select(`[id="rect${key}"]`);
        const newColor = becomeFull ? color(key) : 'white';
        toggle(key, null, newColor);
        thatRect.attr('fill', curData[key].line);
      }
    }
    const fullColor = entry === all ? 'black' : color(entry);
    const newColor = becomeFull ? fullColor : 'white';
    toggle(entry, null, newColor);
    rect.attr('fill', curData[entry].line);
  });
  // dot toggle in the legend
  const dot = legend.append('circle')
    .attr('id', `dot${entry}`)
    .attr('class', 'legend circle')
    .attr('cx', 30)
    .attr('cy', ly)
    .attr('r', 5)
    .attr('stroke', 'darkGrey')
    .attr('fill', curData[entry].dots);
  dot.on('click', () => {
    const becomeFull = curData[entry].dotsOpacity === 0;
    if (entry === all) {
      for (const key of Object.keys(curData).slice(1)) {
        const thatCircle = select(`[id="dot${key}"]`);
        const newColor = becomeFull ? color(key) : 'white';
        toggle(key, newColor, null);
        thatCircle.attr('fill', curData[key].dots);
      }
    }
    const fullColor = entry === all ? 'black' : color(entry);
    const newColor = becomeFull ? fullColor : 'white';
    toggle(entry, newColor, null);
    dot.attr('fill', curData[entry].dots);
  });
  // text for key in the legend
  legend.append('text')
      .attr('class', 'legend')
      .attr('x', 40)
      .attr('y', ly + 5)
      .style('font', '8px sans-serif')
      .text(decodeURI(entry));
}

